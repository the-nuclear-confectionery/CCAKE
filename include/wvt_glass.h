#ifndef WVT_GLASS_H
#define WVT_GLASS_H

#include <string>
#include <vector>

namespace ccake {
/// @brief Adaptive Weighted Voronoi Tessellation (WVT) glass initial conditions.
/// @details Resamples a gridded CCAKE initial condition onto a glass-like particle
/// cloud whose number density tracks the source (Diehl et al., PASA 32 (2015) e048,
/// arXiv:1211.0525). Particles are placed by rejection sampling of the driving density
/// (column 4 of the file), given h_i = eta_smooth (m / rho_i)^(1/D), and relaxed by the
/// repulsive pseudo-force of Diehl Eq. 2. Every source column is then interpolated
/// (bi-/tri-linearly) at the relaxed positions, so each glass row has the same layout as
/// a row of the source file and goes through the normal IC reader.
///
/// The glass is made BEFORE the cutoff: it relaxes over the whole source support (cells
/// with nonzero density), so the energy cutoff surface lies inside the glass rather than
/// on its reflective boundary. CCAKE's readers then apply energy_cutoff to the glass rows
/// exactly as they do to the rows of a grid IC.
///
/// This is the former standalone utilities/wvt_ic_generator (adaptive mode), built into
/// CCAKE: the IC readers consume the glass directly and the 3D stretch defaults to the
/// kernel's eta_h / eta_h_eta.
namespace wvt {

struct GlassSettings
{
  int    number_of_particles = 4000; ///< placed particles, before the cutoff drop
  int    iterations          = 100;  ///< relaxation iterations
  double tolerance           = 0.0;  ///< stop once max|dx|/h < tolerance (0: run all)
  double eta_smooth          = 1.0;  ///< h_i = eta_smooth (m / rho_i)^(1/D); the push reaches 2h.
                                     ///< Larger values drain the core (1.2: core -15%, 1.0: -7%)
  double eta_stretch         = 1.0;  ///< 3D relaxation frame z = eta_stretch * eta
  double mu0                 = 0.1;  ///< initial relaxation gain, annealed to 0
  double softening_fraction  = 0.1;  ///< force softening, in units of h_ij
  double max_move            = 0.5;  ///< per-iteration step cap, in units of h_i
  unsigned seed              = 0;    ///< RNG seed of the initial placement
};

/// The relaxed glass: rows in the source file's column layout plus the grid header
/// the standalone generator wrote for it (mean spacing, eta step, glass minima).
struct Glass
{
  int dimension = 2;
  int number_of_columns = 0;
  std::vector<double> rows;          ///< number_of_rows() * number_of_columns
  double step_x = 0.0, step_y = 0.0, step_eta = 0.0;
  double x_min = 0.0, y_min = 0.0, eta_min = 0.0;

  int number_of_rows() const { return number_of_columns > 0
                                      ? static_cast<int>(rows.size()) / number_of_columns : 0; }
  const double* row(int i) const { return &rows[static_cast<size_t>(i) * number_of_columns]; }
};

/// Reads a gridded CCAKE IC (header "#0 dx dy deta 0 xmin ymin etamin", then rows) and
/// returns its relaxed WVT glass.
Glass make_glass(const std::string& source_file, const GlassSettings& settings);

/// Writes the glass as a CCAKE IC file, byte for byte in the standalone generator's format.
void write_glass(const std::string& path, const Glass& glass);

} // namespace wvt
} // namespace ccake
#endif
