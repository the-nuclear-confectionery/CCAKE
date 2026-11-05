// ============================================================================
// Adaptive WVT glass initial conditions
// (Diehl, Rockefeller, Fryer, Riethmiller, Statler, PASA 32 (2015) e048).
//
// What make_glass() does, in order:
//
//   1. read_source            the gridded IC file -> SourceGrid (all columns per cell,
//                             bi-/tri-linear interpolation between cells)
//   2. find_relaxation_box    the box the particles live in: every cell with nonzero
//                             density plus a margin. NO energy cutoff: the glass is made
//                             before the cutoff, which CCAKE's readers apply afterwards.
//   3. place_particles_by_density
//                             N equal-"mass" particles by rejection sampling of rho
//   4. relax                  repeat `iterations` times:
//                               update_smoothing_lengths   h_i = eta_smooth (m/rho_i)^(1/D)
//                               compute_wvt_displacements  Diehl Eq. 2 (repulsion), with
//                                                          reflective ghosts at the box faces
//                               move_particles             annealed, capped step
//   5. drop_vacuum_particles  particles pushed to density 0 carry nothing
//   6. resample_source_rows   every source column interpolated at the relaxed positions
//   7. set_glass_header       mean spacing and glass minima (the header the standalone
//                             generator wrote; CCAKE sizes its neighbour grid from it)
//
// In 3D the force acts in a stretched frame z = eta_stretch * eta, so the repulsion is
// isotropic there and the glass gets the anisotropy CCAKE's kernel expects.
// ============================================================================
#include "wvt_glass.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>

#include <Kokkos_Core.hpp>
#include <Cabana_Core.hpp>

#include "formatted_output.h"
#include "utilities.h"   // MemorySpace, ExecutionSpace, ListType, sort_neighbor_rows

namespace ccake {
namespace wvt {
namespace {

//==============================================================================
// Basic types and helpers
//==============================================================================

using HostExecution = Kokkos::DefaultHostExecutionSpace;

/// A particle position (x, y, eta). For a 2D source eta is the single slice value.
using Position = std::array<double, 3>;

/// Positions handed to Cabana's neighbour search (stretched frame).
using PositionAoSoA = Cabana::AoSoA<Cabana::MemberTypes<double[3]>, MemorySpace, VECTOR_LENGTH>;

constexpr int MAX_COLUMNS = 32;
constexpr int COLUMN_X = 0, COLUMN_Y = 1, COLUMN_ETA = 2, COLUMN_DENSITY = 3;

constexpr int    RELAXATION_BOX_MARGIN = 4;    ///< grid cells of padding around the support
constexpr int    ENVELOPE_LATTICE      = 80;   ///< lattice points per axis for the sampling envelope
constexpr double ENVELOPE_SAFETY       = 1.05; ///< envelope = ENVELOPE_SAFETY * lattice maximum
constexpr double RHO_FLOOR_FRACTION    = 1e-2; ///< rho is floored at this fraction of the peak
constexpr int    REPORT_EVERY          = 10;   ///< iterations between progress lines

/// A parallel loop over [0, n) on the host (all containers here are host memory).
Kokkos::RangePolicy<HostExecution> host_range(int n)
{
  return Kokkos::RangePolicy<HostExecution>(0, n);
}

[[noreturn]] void abort_with(const std::string& message)
{
  std::cerr << "ERROR (WVT glass): " << message << std::endl;
  Kokkos::finalize();
  abort();
}

//==============================================================================
// The source grid
//==============================================================================

/// The region the particles are confined to (box faces are reflective).
struct RelaxationBox
{
  double x_min = 0, x_max = 0;
  double y_min = 0, y_max = 0;
  double eta_min = 0, eta_max = 0;   ///< both equal the slice value for a 2D source
};

/// The gridded IC: every column at every (ix, iy, ieta) cell, plus the relaxation box.
struct SourceGrid
{
  int number_of_columns = 17;
  int dimension = 2;                  ///< 2: a single eta slice, 3: several
  int nx = 0, ny = 0, neta = 1;       ///< cells per axis
  double x0 = 0, y0 = 0, eta0 = 0;    ///< position of cell (0, 0, 0)
  double dx = 1, dy = 1, deta = 1;    ///< cell size
  double cell_volume = 1;             ///< dx dy deta (3D) or dx dy (2D)
  std::vector<double> values;         ///< nx*ny*neta*number_of_columns, row-major (ix, iy, ieta, column)
  RelaxationBox box;

  bool is_3d() const { return dimension == 3; }

  const double* cell(int ix, int iy, int ieta) const
  {
    return &values[(((static_cast<size_t>(ix) * ny + iy) * neta) + ieta) * number_of_columns];
  }

  /// Tri-linear interpolation of all columns at (x, y, eta); zero outside the grid.
  /// For a 2D source (neta == 1) it is a bi-linear interpolation independent of eta.
  void sample(double x, double y, double eta, double* out) const
  {
    // cell containing the point and the fractional position inside it
    double fx = (x - x0) / dx, fy = (y - y0) / dy;
    int i = static_cast<int>(std::floor(fx));
    int j = static_cast<int>(std::floor(fy));
    if (i < 0 || j < 0 || i >= nx - 1 || j >= ny - 1) {
      for (int c = 0; c < number_of_columns; ++c) out[c] = 0.0;
      return;
    }
    double tx = fx - i, ty = fy - j;

    // eta is clamped to the first/last slice rather than zeroed
    int k = 0; double tz = 0.0;
    if (neta > 1) {
      double fz = (eta - eta0) / deta;
      k = static_cast<int>(std::floor(fz));
      if (k < 0) { k = 0; tz = 0.0; }
      else if (k >= neta - 1) { k = neta - 2; tz = 1.0; }
      else tz = fz - k;
    }
    int k1 = k;
    if (neta > 1) k1 = k + 1;

    const double* c000 = cell(i, j, k);      const double* c100 = cell(i + 1, j, k);
    const double* c010 = cell(i, j + 1, k);  const double* c110 = cell(i + 1, j + 1, k);
    const double* c001 = cell(i, j, k1);     const double* c101 = cell(i + 1, j, k1);
    const double* c011 = cell(i, j + 1, k1); const double* c111 = cell(i + 1, j + 1, k1);
    const double w00 = (1 - tx) * (1 - ty), w10 = tx * (1 - ty),
                 w01 = (1 - tx) * ty,       w11 = tx * ty;
    for (int c = 0; c < number_of_columns; ++c) {
      double lower = w00 * c000[c] + w10 * c100[c] + w01 * c010[c] + w11 * c110[c];
      double upper = w00 * c001[c] + w10 * c101[c] + w01 * c011[c] + w11 * c111[c];
      out[c] = (1 - tz) * lower + tz * upper;
    }
  }

  /// The driving density (column 4) at (x, y, eta).
  double density(double x, double y, double eta) const
  {
    double out[MAX_COLUMNS];
    sample(x, y, eta, out);
    return out[COLUMN_DENSITY];
  }

  /// Integral of the driving density over the grid. Deliberately serial: it sets the
  /// particle mass and hence every h_i, so a reduction in unspecified order would make
  /// the whole glass depend on the thread count.
  double total_mass() const
  {
    double sum = 0.0;
    for (int ix = 0; ix < nx; ++ix)
      for (int iy = 0; iy < ny; ++iy)
        for (int ieta = 0; ieta < neta; ++ieta) sum += cell(ix, iy, ieta)[COLUMN_DENSITY];
    return sum * cell_volume;
  }

  /// Largest driving density of any cell.
  double peak_density() const
  {
    double peak = 0.0;
    for (int ix = 0; ix < nx; ++ix)
      for (int iy = 0; iy < ny; ++iy)
        for (int ieta = 0; ieta < neta; ++ieta) peak = std::max(peak, cell(ix, iy, ieta)[COLUMN_DENSITY]);
    return peak;
  }
};

//==============================================================================
// Step 1: reading the source file
//==============================================================================

/// One line of the file buffer, [begin, end) without the newline.
struct LineSpan
{
  const char* begin;
  const char* end;
};

std::string read_whole_file(const std::string& path)
{
  std::ifstream in(path, std::ios::binary);
  if (!in) abort_with("cannot open source file " + path);
  in.seekg(0, std::ios::end);
  const std::streamoff file_size = in.tellg();
  in.seekg(0, std::ios::beg);
  std::string buffer(static_cast<size_t>(file_size), '\0');
  if (file_size > 0) in.read(&buffer[0], file_size);
  return buffer;
}

/// Returns the line starting at `cursor` and advances `cursor` past its newline.
LineSpan next_line(const char*& cursor, const char* buffer_end)
{
  LineSpan line;
  line.begin = cursor;
  while (cursor < buffer_end && *cursor != '\n') ++cursor;
  line.end = cursor;
  if (cursor < buffer_end) ++cursor;
  return line;
}

/// Reads the grid header "#0 dx dy deta 0 xmin ymin etamin". The eta origin is taken
/// from the data instead (set_grid_layout), so the header's etamin is not used.
void parse_grid_header(const LineSpan& header, SourceGrid& source)
{
  std::string line(header.begin, header.end);
  size_t hash = line.find('#');
  if (hash != std::string::npos) line = line.substr(hash + 1);
  std::istringstream header_stream(line);
  double ignore, eta_min_of_header = 0.0;
  header_stream >> ignore >> source.dx >> source.dy >> source.deta >> ignore
                >> source.x0 >> source.y0 >> eta_min_of_header;
}

/// The data lines after the header: blank lines and '#' comment lines are skipped.
std::vector<LineSpan> collect_data_lines(const char* cursor, const char* buffer_end,
                                         size_t expected_lines)
{
  std::vector<LineSpan> lines;
  lines.reserve(expected_lines);
  while (cursor < buffer_end) {
    LineSpan line = next_line(cursor, buffer_end);
    if (line.begin == line.end || *line.begin == '#') continue;
    lines.push_back(line);
  }
  return lines;
}

/// Parses each line into up to MAX_COLUMNS numbers, in parallel (a strtod walk is far
/// cheaper than one istringstream per line on a ~150 MB source).
void parse_rows(const std::vector<LineSpan>& lines,
                std::vector<std::array<double, MAX_COLUMNS>>& rows,
                std::vector<int>& columns_in_row)
{
  const int number_of_rows = static_cast<int>(lines.size());
  rows.resize(number_of_rows);
  columns_in_row.resize(number_of_rows);
  Kokkos::parallel_for("wvt_parse_rows", host_range(number_of_rows), [&](const int row) {
    const char* p = lines[row].begin;
    const char* line_end = lines[row].end;
    int n = 0;
    while (n < MAX_COLUMNS && p < line_end) {
      char* next = nullptr;
      double value = std::strtod(p, &next);
      if (next == p) break;
      rows[row][n++] = value;
      p = next;
    }
    columns_in_row[row] = n;
  });
}

/// Sets the dimension, the eta origin, the cell counts and the cell volume from the
/// data, and returns each row's integer cell index (ix, iy, ieta).
std::vector<std::array<int, 3>> set_grid_layout(
    const std::vector<std::array<double, MAX_COLUMNS>>& rows, SourceGrid& source)
{
  const int number_of_rows = static_cast<int>(rows.size());

  // Dimensionality from the DATA eta range: a 2D IC is a single eta slice, which the
  // header may place at a nonzero index. The eta grid starts at the data minimum.
  double eta_lowest = 1e30, eta_highest = -1e30;
  for (int row = 0; row < number_of_rows; ++row) {
    eta_lowest = std::min(eta_lowest, rows[row][COLUMN_ETA]);
    eta_highest = std::max(eta_highest, rows[row][COLUMN_ETA]);
  }
  const bool is_3d = (source.deta > 0.0) && (eta_highest - eta_lowest) > 0.5 * source.deta;
  source.eta0 = eta_lowest;

  int ix_max = 0, iy_max = 0, ieta_max = 0;
  std::vector<std::array<int, 3>> cell_index(number_of_rows);
  for (int row = 0; row < number_of_rows; ++row) {
    int ix = static_cast<int>(std::lround((rows[row][COLUMN_X] - source.x0) / source.dx));
    int iy = static_cast<int>(std::lround((rows[row][COLUMN_Y] - source.y0) / source.dy));
    int ieta = 0;
    if (is_3d) ieta = static_cast<int>(std::lround((rows[row][COLUMN_ETA] - source.eta0) / source.deta));
    cell_index[row] = {ix, iy, ieta};
    ix_max = std::max(ix_max, ix);
    iy_max = std::max(iy_max, iy);
    ieta_max = std::max(ieta_max, ieta);
  }

  source.nx = ix_max + 1;
  source.ny = iy_max + 1;
  if (is_3d) {
    source.neta = ieta_max + 1;
    source.dimension = 3;
    source.cell_volume = source.dx * source.dy * source.deta;
  } else {
    source.neta = 1;
    source.dimension = 2;
    source.cell_volume = source.dx * source.dy;
  }
  return cell_index;
}

/// Copies every row into its cell.
void fill_grid(const std::vector<std::array<double, MAX_COLUMNS>>& rows,
               const std::vector<std::array<int, 3>>& cell_index, SourceGrid& source)
{
  const int number_of_rows = static_cast<int>(rows.size());
  source.values.assign(static_cast<size_t>(source.nx) * source.ny * source.neta
                       * source.number_of_columns, 0.0);
  for (int row = 0; row < number_of_rows; ++row) {
    const size_t flat_cell = (static_cast<size_t>(cell_index[row][0]) * source.ny
                              + cell_index[row][1]) * source.neta + cell_index[row][2];
    double* destination = &source.values[flat_cell * source.number_of_columns];
    for (int c = 0; c < source.number_of_columns; ++c) destination[c] = rows[row][c];
  }
}

/// Reads a gridded CCAKE IC file into a SourceGrid (the box is set separately).
SourceGrid read_source(const std::string& path)
{
  const std::string buffer = read_whole_file(path);
  const char* const buffer_end = buffer.data() + buffer.size();
  const char* cursor = buffer.data();

  SourceGrid source;
  parse_grid_header(next_line(cursor, buffer_end), source);

  const std::vector<LineSpan> lines = collect_data_lines(cursor, buffer_end, buffer.size() / 40 + 16);
  if (lines.empty()) abort_with("no data rows in " + path);

  std::vector<std::array<double, MAX_COLUMNS>> rows;
  std::vector<int> columns_in_row;
  parse_rows(lines, rows, columns_in_row);
  source.number_of_columns = columns_in_row[0];

  const std::vector<std::array<int, 3>> cell_index = set_grid_layout(rows, source);
  fill_grid(rows, cell_index, source);
  return source;
}

//==============================================================================
// Step 2: the relaxation box
//==============================================================================

/// Every cell with nonzero density (the whole support, no energy cutoff), padded by
/// RELAXATION_BOX_MARGIN cells and clipped to the grid.
RelaxationBox find_relaxation_box(const SourceGrid& source)
{
  int ix_first = source.nx, ix_last = 0;
  int iy_first = source.ny, iy_last = 0;
  int ieta_first = source.neta, ieta_last = 0;
  for (int ix = 0; ix < source.nx; ++ix)
    for (int iy = 0; iy < source.ny; ++iy)
      for (int ieta = 0; ieta < source.neta; ++ieta)
        if (source.cell(ix, iy, ieta)[COLUMN_DENSITY] > 0.0) {
          ix_first = std::min(ix_first, ix);       ix_last = std::max(ix_last, ix);
          iy_first = std::min(iy_first, iy);       iy_last = std::max(iy_last, iy);
          ieta_first = std::min(ieta_first, ieta); ieta_last = std::max(ieta_last, ieta);
        }
  if (ix_first > ix_last) {   // no positive cell at all: use the whole grid
    ix_first = 0;   ix_last = source.nx - 1;
    iy_first = 0;   iy_last = source.ny - 1;
    ieta_first = 0; ieta_last = source.neta - 1;
  }

  const int margin = RELAXATION_BOX_MARGIN;
  ix_first = std::max(ix_first - margin, 0);     ix_last = std::min(ix_last + margin, source.nx - 1);
  iy_first = std::max(iy_first - margin, 0);     iy_last = std::min(iy_last + margin, source.ny - 1);
  ieta_first = std::max(ieta_first - margin, 0); ieta_last = std::min(ieta_last + margin, source.neta - 1);

  RelaxationBox box;
  box.x_min = source.x0 + ix_first * source.dx;     box.x_max = source.x0 + ix_last * source.dx;
  box.y_min = source.y0 + iy_first * source.dy;     box.y_max = source.y0 + iy_last * source.dy;
  box.eta_min = source.eta0 + ieta_first * source.deta;
  box.eta_max = source.eta0 + ieta_last * source.deta;
  if (!source.is_3d()) { box.eta_min = source.eta0; box.eta_max = source.eta0; }
  return box;
}

//==============================================================================
// Step 3: initial placement
//==============================================================================

/// Upper bound on the density inside the box for rejection sampling: the maximum on an
/// (ENVELOPE_LATTICE+1)^D lattice times ENVELOPE_SAFETY. (An exact max, so the
/// parallel order does not matter.)
double rejection_envelope(const SourceGrid& source)
{
  const RelaxationBox& box = source.box;
  const int M = ENVELOPE_LATTICE;
  int eta_points = 0;
  if (source.is_3d()) eta_points = M;
  double density_max = 0.0;
  Kokkos::parallel_reduce("wvt_peak_scan", host_range((M + 1) * (M + 1)),
    [&](const int ab, double& local_max) {
      const int a = ab / (M + 1), b = ab % (M + 1);
      for (int c = 0; c <= eta_points; ++c) {
        const double x = box.x_min + (box.x_max - box.x_min) * a / M;
        const double y = box.y_min + (box.y_max - box.y_min) * b / M;
        double eta = source.eta0;
        if (source.is_3d()) eta = box.eta_min + (box.eta_max - box.eta_min) * c / std::max(M, 1);
        local_max = std::max(local_max, source.density(x, y, eta));
      }
    }, Kokkos::Max<double>(density_max));
  density_max *= ENVELOPE_SAFETY;
  return density_max;
}

/// N positions distributed proportionally to the density inside the box (rejection
/// sampling; serial, so the result depends only on the seed).
std::vector<Position> place_particles_by_density(int number_of_particles, const SourceGrid& source,
                                                 std::mt19937_64& rng)
{
  const RelaxationBox& box = source.box;
  const double envelope = rejection_envelope(source);
  std::uniform_real_distribution<double> random_x(box.x_min, box.x_max);
  std::uniform_real_distribution<double> random_y(box.y_min, box.y_max);
  std::uniform_real_distribution<double> random_eta(box.eta_min, box.eta_max);
  std::uniform_real_distribution<double> random_acceptance(0.0, envelope);

  std::vector<Position> positions;
  positions.reserve(number_of_particles);
  while (static_cast<int>(positions.size()) < number_of_particles) {
    // the draw order (x, y, [eta], acceptance) fixes the glass for a given seed
    double x = random_x(rng);
    double y = random_y(rng);
    double eta = source.eta0;
    if (source.is_3d()) eta = random_eta(rng);
    if (source.density(x, y, eta) > random_acceptance(rng)) positions.push_back({x, y, eta});
  }
  return positions;
}

//==============================================================================
// Step 4: the relaxation
//==============================================================================

/// h_i = eta_smooth (m / rho_i)^(1/D): each particle "carries" the same integrated density
/// m, so the number density tracks rho. rho is floored at rho_floor, which caps h.
void update_smoothing_lengths(const std::vector<Position>& positions, const SourceGrid& source,
                              double mass_per_particle, double eta_smooth, double rho_floor,
                              std::vector<double>& h)
{
  const int number_of_particles = static_cast<int>(positions.size());
  h.resize(number_of_particles);
  const double inverse_dimension = 1.0 / source.dimension;
  Kokkos::parallel_for("wvt_update_h", host_range(number_of_particles), [&](const int i) {
    double rho = std::max(source.density(positions[i][0], positions[i][1], positions[i][2]), rho_floor);
    h[i] = eta_smooth * std::pow(mass_per_particle / rho, inverse_dimension);
  });
}

/// Mirror images of the particles near the box faces. They only push (never move), so
/// edge particles are not drawn outward by missing neighbours. Stored stretched.
struct GhostParticles
{
  std::vector<Position> stretched_positions;
  std::vector<double> h;
};

/// Where the ghosts are made: the box faces and the band around them.
struct GhostFaces
{
  RelaxationBox box;
  double band;            ///< mirror particles within this distance of a face (2 h_max)
  double z_min, z_max;    ///< the eta faces in the stretched frame
  double eta_stretch;
  bool is_3d;
};

/// Appends the ghosts of particle p, in the fixed face order x_min, x_max, y_min, y_max,
/// eta_min, eta_max.
void add_ghosts_of(const Position& p, double h, const GhostFaces& faces, GhostParticles& ghosts)
{
  const double x = p[0], y = p[1], z = faces.eta_stretch * p[2];
  const RelaxationBox& box = faces.box;
  Position mirrors[6];
  int count = 0;
  if (x < box.x_min + faces.band) mirrors[count++] = {2 * box.x_min - x, y, z};
  if (x > box.x_max - faces.band) mirrors[count++] = {2 * box.x_max - x, y, z};
  if (y < box.y_min + faces.band) mirrors[count++] = {x, 2 * box.y_min - y, z};
  if (y > box.y_max - faces.band) mirrors[count++] = {x, 2 * box.y_max - y, z};
  if (faces.is_3d && z < faces.z_min + faces.band) mirrors[count++] = {x, y, 2 * faces.z_min - z};
  if (faces.is_3d && z > faces.z_max - faces.band) mirrors[count++] = {x, y, 2 * faces.z_max - z};
  for (int m = 0; m < count; ++m) {
    ghosts.stretched_positions.push_back(mirrors[m]);
    ghosts.h.push_back(h);
  }
}

/// All ghosts, in particle order. The order matters: it is the order in which ghosts
/// enter the neighbour lists and hence the force sums. A serial loop gives it directly,
/// and it is cheap next to the neighbour search.
GhostParticles make_reflective_ghosts(const std::vector<Position>& positions,
                                      const std::vector<double>& h, const GhostFaces& faces)
{
  GhostParticles ghosts;
  const int number_of_particles = static_cast<int>(positions.size());
  for (int i = 0; i < number_of_particles; ++i) add_ghosts_of(positions[i], h[i], faces, ghosts);
  return ghosts;
}

/// Bounds of Cabana's background binning grid: the box plus 2.5 h_max, snapped to whole
/// cells of size 2 h_max with a one-cell pad. Non-round bounds (e.g. an irrational
/// stretch) make Cabana's bin count degenerate.
void neighbor_grid_bounds(const GhostFaces& faces, double h_max,
                          double min_position[3], double max_position[3])
{
  const RelaxationBox& box = faces.box;
  min_position[0] = box.x_min - 2.5 * h_max; max_position[0] = box.x_max + 2.5 * h_max;
  min_position[1] = box.y_min - 2.5 * h_max; max_position[1] = box.y_max + 2.5 * h_max;
  if (faces.is_3d) {
    min_position[2] = faces.z_min - 2.5 * h_max;
    max_position[2] = faces.z_max + 2.5 * h_max;
  } else {
    min_position[2] = -h_max;
    max_position[2] =  h_max;
  }
  const double cell = 2.0 * h_max;
  for (int idir = 0; idir < 3; ++idir) {
    min_position[idir] = std::floor(min_position[idir] / cell - 1.0) * cell;
    max_position[idir] = std::ceil (max_position[idir] / cell + 1.0) * cell;
  }
}

/// The repulsive pseudo-force on every real particle i (ghosts only push):
///   F_i = sum_j h_ij f(h_ij, r_ij) rhat_ij,   h_ij = (h_i + h_j)/2,
///   f(h, r) = (h/(r+eps))^2 - (h/(2h+eps))^2  for r < 2h,  eps = softening_fraction h.
/// f is Diehl's Eq. 2; the prefactor is the symmetrised h_ij of Arth et al. 2019
/// (arXiv:1907.11250, Sec. 3.1) instead of Diehl's h_i, so i and j push each other equally.
/// Diehl's h_i lets a dilute-halo particle (large h) push harder than a core particle
/// pushes back, which drains the core. The second term makes f vanish at r = 2h.
/// Positions are in the stretched frame.
Kokkos::View<double*[3], MemorySpace> repulsive_forces(
    const PositionAoSoA& stretched, const Kokkos::View<double*, MemorySpace>& h_view,
    const ListType& list, int number_of_particles, double softening_fraction, bool is_3d)
{
  auto position = Cabana::slice<0>(stretched);
  Kokkos::View<double*[3], MemorySpace> force("wvt_force", number_of_particles);
  Kokkos::parallel_for(
      "wvt_force", Kokkos::RangePolicy<ExecutionSpace>(0, number_of_particles),
      KOKKOS_LAMBDA(const int i) {
        double xi = position(i, 0), yi = position(i, 1), zi = position(i, 2), hi = h_view(i);
        const int number_of_neighbors = Cabana::NeighborList<ListType>::numNeighbor(list, i);
        double fx = 0.0, fy = 0.0, fz = 0.0;
        for (int n = 0; n < number_of_neighbors; ++n) {
          int j = Cabana::NeighborList<ListType>::getNeighbor(list, i, n);
          double dxx = xi - position(j, 0), dyy = yi - position(j, 1);
          double dzz = 0.0;
          if (is_3d) dzz = zi - position(j, 2);
          double r = std::sqrt(dxx * dxx + dyy * dyy + dzz * dzz);
          if (r <= 0.0) continue;
          double hij = 0.5 * (hi + h_view(j));
          if (r >= 2.0 * hij) continue;
          double eps = softening_fraction * hij;
          double a = hij / (r + eps);
          double b = hij / (2.0 * hij + eps);
          double f = a * a - b * b;
          fx += hij * f * dxx / r; fy += hij * f * dyy / r; fz += hij * f * dzz / r;
        }
        force(i, 0) = fx; force(i, 1) = fy; force(i, 2) = fz;
      });
  Kokkos::fence();
  return force;
}

/// One WVT displacement per particle (before the gain mu is applied): reflective ghosts,
/// a per-particle-radius neighbour list over reals + ghosts, the repulsive force, and the
/// eta component un-stretched back to rapidity.
void compute_wvt_displacements(const std::vector<Position>& positions, const std::vector<double>& h,
                               const SourceGrid& source, double softening_fraction,
                               double eta_stretch, std::vector<Position>& displacement)
{
  const int number_of_particles = static_cast<int>(positions.size());
  double h_max = 0.0;
  Kokkos::parallel_reduce("wvt_hmax", host_range(number_of_particles),
    [&](const int i, double& local_max) { local_max = std::max(local_max, h[i]); },
    Kokkos::Max<double>(h_max));

  GhostFaces faces;
  faces.box = source.box;
  faces.band = 2.0 * h_max;
  faces.z_min = eta_stretch * source.box.eta_min;
  faces.z_max = eta_stretch * source.box.eta_max;
  faces.eta_stretch = eta_stretch;
  faces.is_3d = source.is_3d();
  const GhostParticles ghosts = make_reflective_ghosts(positions, h, faces);
  const int number_of_ghosts = static_cast<int>(ghosts.h.size());
  const int number_with_ghosts = number_of_particles + number_of_ghosts;

  // reals first, then ghosts; positions stretched, search radius 2 h per particle
  PositionAoSoA stretched("wvt_positions", number_with_ghosts);
  auto position = Cabana::slice<0>(stretched);
  Kokkos::View<double*, MemorySpace> h_view("wvt_h", number_with_ghosts);
  Kokkos::View<double*, MemorySpace> search_radius("wvt_radii", number_with_ghosts);
  Kokkos::parallel_for("wvt_fill_real", host_range(number_of_particles), [&](const int i) {
    position(i, 0) = positions[i][0];
    position(i, 1) = positions[i][1];
    position(i, 2) = eta_stretch * positions[i][2];
    h_view(i) = h[i];
    search_radius(i) = 2.0 * h[i];
  });
  Kokkos::parallel_for("wvt_fill_ghost", host_range(number_of_ghosts), [&](const int g) {
    const int i = number_of_particles + g;
    position(i, 0) = ghosts.stretched_positions[g][0];
    position(i, 1) = ghosts.stretched_positions[g][1];
    position(i, 2) = ghosts.stretched_positions[g][2];
    h_view(i) = ghosts.h[g];
    search_radius(i) = 2.0 * ghosts.h[g];
  });

  // Per-particle radii (the list type of CCAKE's per_particle backend): each particle
  // searches its own 2 h_i; 2 h_max only sizes the background bins.
  double min_position[3], max_position[3];
  neighbor_grid_bounds(faces, h_max, min_position, max_position);
  ListType list(position, 0, number_with_ghosts, 2.0 * h_max, search_radius, 1.0,
                min_position, max_position);
  sort_neighbor_rows(list, number_of_particles);   // canonical order of every force sum

  Kokkos::View<double*[3], MemorySpace> force =
      repulsive_forces(stretched, h_view, list, number_of_particles, softening_fraction, faces.is_3d);

  displacement.resize(number_of_particles);
  const bool is_3d = faces.is_3d;
  Kokkos::parallel_for("wvt_displacement_copy", host_range(number_of_particles), [&](const int i) {
    double eta_displacement = 0.0;
    if (is_3d) eta_displacement = force(i, 2) / eta_stretch;   // back to rapidity
    displacement[i] = {force(i, 0), force(i, 1), eta_displacement};
  });
}

/// Applies x_i += mu * displacement_i, with the step capped at max_move * h_i (eta
/// measured in the stretched frame) and the particle clamped into the box. Returns the
/// largest step in units of h, the convergence monitor.
double move_particles(std::vector<Position>& positions, const std::vector<Position>& displacement,
                      const std::vector<double>& h, double mu, const GlassSettings& settings,
                      const SourceGrid& source)
{
  const int number_of_particles = static_cast<int>(positions.size());
  const RelaxationBox& box = source.box;
  const bool is_3d = source.is_3d();
  const double eta_stretch = settings.eta_stretch;
  double largest_step_over_h = 0.0;
  Kokkos::parallel_reduce("wvt_move", host_range(number_of_particles),
    [&](const int i, double& local_largest) {
      double sx = mu * displacement[i][0], sy = mu * displacement[i][1], se = mu * displacement[i][2];
      double step = std::sqrt(sx * sx + sy * sy + (eta_stretch * se) * (eta_stretch * se));
      double cap = settings.max_move * h[i];
      if (step > cap && step > 0.0) { double f = cap / step; sx *= f; sy *= f; se *= f; step = cap; }
      positions[i][0] = std::min(std::max(positions[i][0] + sx, box.x_min), box.x_max);
      positions[i][1] = std::min(std::max(positions[i][1] + sy, box.y_min), box.y_max);
      if (is_3d) positions[i][2] = std::min(std::max(positions[i][2] + se, box.eta_min), box.eta_max);
      local_largest = std::max(local_largest, step / std::max(h[i], 1e-30));
    }, Kokkos::Max<double>(largest_step_over_h));
  return largest_step_over_h;
}

void report_iteration(int iteration, double mu, double largest_step_over_h)
{
  std::ostringstream line;
  line << "WVT iteration " << iteration << "  mu=" << mu << "  max|dx|/h=" << largest_step_over_h;
  formatted_output::detail(line.str());
}

/// The annealed relaxation: the gain mu falls linearly from mu0 to 0 over `iterations`.
/// Stops early once the largest step drops below `tolerance` h (if tolerance > 0); the
/// schedule is not stretched, so an early stop happens mid-anneal.
void relax(std::vector<Position>& positions, const SourceGrid& source,
           const GlassSettings& settings, double mass_per_particle, double rho_floor)
{
  const int iterations = settings.iterations;
  std::vector<double> h;
  std::vector<Position> displacement;
  for (int k = 0; k < iterations; ++k) {
    const double mu = settings.mu0 * (1.0 - static_cast<double>(k) / iterations);
    update_smoothing_lengths(positions, source, mass_per_particle, settings.eta_smooth, rho_floor, h);
    compute_wvt_displacements(positions, h, source, settings.softening_fraction,
                              settings.eta_stretch, displacement);
    const double largest_step_over_h = move_particles(positions, displacement, h, mu, settings, source);

    const bool converged = (settings.tolerance > 0.0) && (largest_step_over_h < settings.tolerance);
    if (k % REPORT_EVERY == 0 || k == iterations - 1 || converged)
      report_iteration(k, mu, largest_step_over_h);
    if (converged) {
      formatted_output::detail("WVT converged: max|dx|/h < " + std::to_string(settings.tolerance)
                               + " at iteration " + std::to_string(k));
      break;
    }
  }
}

//==============================================================================
// Steps 5-7: from relaxed positions to glass rows
//==============================================================================

/// Drops particles at exactly zero density (e.g. pushed into an empty margin cell): their
/// rows would carry nothing. The energy cutoff is NOT applied here. Order is preserved.
std::vector<Position> drop_vacuum_particles(const std::vector<Position>& positions,
                                            const SourceGrid& source)
{
  std::vector<Position> kept;
  kept.reserve(positions.size());
  for (const Position& p : positions)
    if (source.density(p[0], p[1], p[2]) > 0.0) kept.push_back(p);
  return kept;
}

/// One row per particle in the source's column layout: every column interpolated at the
/// particle, then x, y, eta overwritten with its exact position (eta = 0 in 2D).
void resample_source_rows(const std::vector<Position>& positions, const SourceGrid& source,
                          Glass& glass)
{
  const int number_of_particles = static_cast<int>(positions.size());
  const int number_of_columns = source.number_of_columns;
  const bool is_3d = source.is_3d();
  glass.dimension = source.dimension;
  glass.number_of_columns = number_of_columns;
  glass.rows.assign(static_cast<size_t>(number_of_particles) * number_of_columns, 0.0);
  for (int i = 0; i < number_of_particles; ++i) {
    double* row = &glass.rows[static_cast<size_t>(i) * number_of_columns];
    double sampled[MAX_COLUMNS];
    source.sample(positions[i][0], positions[i][1], positions[i][2], sampled);
    for (int c = 0; c < number_of_columns; ++c) row[c] = sampled[c];
    row[COLUMN_X] = positions[i][0];
    row[COLUMN_Y] = positions[i][1];
    row[COLUMN_ETA] = 0.0;
    if (is_3d) row[COLUMN_ETA] = positions[i][2];
  }
}

/// The header the standalone generator wrote: dx = dy = the mean spacing, deta of the
/// source, and the glass's own minima (the eta slice value in 2D).
void set_glass_header(const std::vector<Position>& positions, double mean_spacing,
                      const SourceGrid& source, Glass& glass)
{
  double x_min = 1e30, y_min = 1e30, eta_min = 1e30;
  for (const Position& p : positions) {
    x_min = std::min(x_min, p[0]);
    y_min = std::min(y_min, p[1]);
    eta_min = std::min(eta_min, p[2]);
  }
  if (!source.is_3d()) eta_min = source.eta0;
  glass.step_x = mean_spacing;
  glass.step_y = mean_spacing;
  glass.step_eta = source.deta;
  glass.x_min = x_min;
  glass.y_min = y_min;
  glass.eta_min = eta_min;
}

/// Mean particle spacing (V/N)^(1/D) over the relaxation box.
double mean_spacing_of(const SourceGrid& source, int number_of_particles)
{
  const RelaxationBox& box = source.box;
  const int n = std::max(number_of_particles, 1);
  if (source.is_3d()) {
    const double volume = (box.x_max - box.x_min) * (box.y_max - box.y_min) * (box.eta_max - box.eta_min);
    return std::cbrt(volume / n);
  }
  const double area = (box.x_max - box.x_min) * (box.y_max - box.y_min) * 1.0;
  return std::sqrt(area / n);
}

void report_settings(const std::string& source_file, const SourceGrid& source,
                     const GlassSettings& settings, double peak)
{
  std::ostringstream line;
  line << "WVT glass from " << source_file;
  formatted_output::update(line.str());
  line.str("");
  line << source.dimension << "D, " << source.number_of_columns << " columns, "
       << settings.number_of_particles << " particles, " << settings.iterations << " iterations";
  formatted_output::detail(line.str());
  line.str("");
  line << "eta_stretch " << settings.eta_stretch << ", peak density " << peak
       << " (no cutoff: the readers apply energy_cutoff to the glass)";
  formatted_output::detail(line.str());
}

} // namespace

//==============================================================================
// Public interface
//==============================================================================

Glass make_glass(const std::string& source_file, const GlassSettings& settings)
{
  // 1-2: the source and the region the glass fills
  SourceGrid source = read_source(source_file);
  source.box = find_relaxation_box(source);
  const double peak = source.peak_density();
  report_settings(source_file, source, settings, peak);

  // equal-mass particles; rho_floor caps h. The floor must come from the TRUE peak: an
  // underestimate inflates h_max until the 2 h_max background cell covers the core.
  const int number_of_particles = settings.number_of_particles;
  const double mass_per_particle = source.total_mass() / std::max(number_of_particles, 1);
  const double rho_floor = RHO_FLOOR_FRACTION * peak;

  // 3-4: place and relax
  std::mt19937_64 rng(settings.seed);
  std::vector<Position> positions = place_particles_by_density(number_of_particles, source, rng);
  relax(positions, source, settings, mass_per_particle, rho_floor);

  // 5-7: the glass rows and their header
  const std::vector<Position> kept = drop_vacuum_particles(positions, source);
  Glass glass;
  resample_source_rows(kept, source, glass);
  set_glass_header(kept, mean_spacing_of(source, number_of_particles), source, glass);

  formatted_output::update("WVT glass: " + std::to_string(glass.number_of_rows()) + " of "
                           + std::to_string(number_of_particles)
                           + " particles (vacuum dropped), cutoff follows in the reader");
  return glass;
}

/// Same format as the standalone generator: header at precision 8 in scientific notation,
/// then every column as %.8e, space separated.
void write_glass(const std::string& path, const Glass& glass)
{
  std::ofstream out(path);
  if (!out) abort_with("cannot write glass file " + path);
  out.setf(std::ios::scientific);
  out.precision(8);
  out << "#0 " << glass.step_x << " " << glass.step_y << " " << glass.step_eta << " 0 "
      << glass.x_min << " " << glass.y_min << " " << glass.eta_min << "\n";

  std::string line;
  char field[32];
  for (int i = 0; i < glass.number_of_rows(); ++i) {
    line.clear();
    const double* row = glass.row(i);
    for (int c = 0; c < glass.number_of_columns; ++c) {
      const int length = std::snprintf(field, sizeof(field), "%.8e", row[c]);
      line.append(field, length);
      if (c + 1 < glass.number_of_columns) line += ' ';
      else line += '\n';
    }
    out << line;
  }
  formatted_output::update("WVT glass written to " + path);
}

} // namespace wvt
} // namespace ccake
