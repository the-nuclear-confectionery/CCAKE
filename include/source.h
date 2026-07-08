#ifndef SOURCE_H
#define SOURCE_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <Cabana_Core.hpp>

#include "constants.h"
#include "eom_default.h"
#include "eom_cartesian.h"
#include "formatted_output.h"
#include "jets.h"
#include "kernel.h"
#include "settings.h"
#include "system_state.h"
#include "thermodynamic_info.h"



namespace ccake
{

#ifndef CONCAT
  #define CONCAT_(a, b) a##b
  #define CONCAT(a, b) CONCAT_(a, b)
#endif
  using SRC = Cabana::MemberTypes<double[3], // position
                                  double[3], // momentum
                                  double,    // energy
                                  double,    // total energy
                                  double,    // mass
                                  double,    // baryon charge
                                  double,    // strangeness charge
                                  double,    // electric charge
                                  double,    // time
                                  int>;       // PID

  namespace SRC_enum{
    enum SRC_members
    {
      position = 0,
      momentum,
      energy,
      total_energy,
      mass,
      baryon_charge,
      strangeness_charge,
      electric_charge,
      time,
      PID
    };
  }

  using source_array = Cabana::AoSoA<SRC, DeviceType,VECTOR_LENGTH>;

  #define SRC_VIEW(prefix, SRC_aosoa) \
  auto CONCAT(prefix, position) = Cabana::slice<SRC_enum::position>(SRC_aosoa); \
  auto CONCAT(prefix, momentum) = Cabana::slice<SRC_enum::momentum>(SRC_aosoa); \
  auto CONCAT(prefix, energy) = Cabana::slice<SRC_enum::energy>(SRC_aosoa); \
  auto CONCAT(prefix, total_energy) = Cabana::slice<SRC_enum::total_energy>(SRC_aosoa); \
  auto CONCAT(prefix, mass) = Cabana::slice<SRC_enum::mass>(SRC_aosoa); \
  auto CONCAT(prefix, baryon_charge) = Cabana::slice<SRC_enum::baryon_charge>(SRC_aosoa); \
  auto CONCAT(prefix, strangeness_charge) = Cabana::slice<SRC_enum::strangeness_charge>(SRC_aosoa); \
  auto CONCAT(prefix, electric_charge) = Cabana::slice<SRC_enum::electric_charge>(SRC_aosoa); \
  auto CONCAT(prefix, time) = Cabana::slice<SRC_enum::time>(SRC_aosoa); \
  auto CONCAT(prefix, PID) = Cabana::slice<SRC_enum::PID>(SRC_aosoa); \



/// @brief External source terms (energy/momentum/charge currents) deposited into
/// the fluid. Handles IC sources read from file (AMPT partons, SMASH particles),
/// the JEWEL scattering-centre reader, and the real-time jet sources (toy model,
/// BBMG). Owned and driven by SPHWorkstation.
template <unsigned int D>
class Source
{
private:

  bool has_jet_sources = false;         ///< real-time jet source (toy model / BBMG)
  bool has_propagating_sources = false; ///< propagating/depleting IC sources (SMASH/AMPT)

  // Persistent master state for propagating IC sources. The launch data is
  // immutable; src_remaining is the per-source energy budget (depleted each
  // step) and src_charge_dumped tracks whether the one-shot charge release
  // has already happened on absorption. The per-step `sources` AoSoA is rebuilt
  // from this master each step (mirrors the BBMG jetInfo -> sources split).
  // The src_* names are generic: they serve BOTH SMASH particles and AMPT partons.
  std::vector<double> src_x, src_y, src_z;          ///< current (propagating) position
  std::vector<double> src_px, src_py, src_pz;       ///< momentum (immutable, Cartesian)
  std::vector<double> src_mass, src_E0;             ///< mass and launch energy (immutable)
  std::vector<double> src_baryon, src_strange, src_charge; ///< charges to dump on absorption
  std::vector<double> src_formation_time;           ///< activation time t_src (immutable)
  std::vector<double> src_remaining;                ///< remaining energy budget (mutable)
  std::vector<char>   src_charge_dumped;            ///< whether charge already released
  std::vector<int>    src_PID;                      ///< particle id

  // Per-timestep deposit state, computed ONCE per timestep in advance_sources()
  // and only READ by build_step_sources() during the RK substeps. This prevents
  // the source from advancing/depleting multiple times per timestep
  // (get_time_derivatives runs once per RK substep).
  std::vector<double> src_dE;                        ///< energy chunk to deposit this step
  std::vector<double> src_vx, src_vy, src_veta;      ///< Milne velocity at deposit (for momentum)
  std::vector<char>   src_dump_now;                  ///< release full charge this step (absorption)

  // Diagnostics: lifetime of each propagating source (steps from activation to
  // absorption) and whether any source escapes with energy/charge left.
  std::vector<int>    src_activate_step;             ///< step index when first active (-1 = never)
  std::vector<int>    src_absorb_step;               ///< step index when absorbed   (-1 = never)
  int                 src_step = 0;                   ///< advance_sources call counter

  double total_energy_transferred;             ///< Total energy transferred by the source term
  double total_baryon_charge_transferred;      ///< Total baryon charge transferred
  double total_strangeness_charge_transferred; ///< Total strangeness charge transferred
  double total_electric_charge_transferred;    ///< Total electric charge transferred

  std::size_t n_sources = 0; ///< Number of sources in the simulation

  source_array sources; ///< Array of source terms for the particles

  std::shared_ptr<Settings> settingsPtr    = nullptr; ///< Settings parsed from input file
  std::shared_ptr<SystemState<D>> systemPtr = nullptr; ///< SPH system (linked list, particles, ...)

  BBMG<D>* bbmgPtr = nullptr; ///< Non-owning ref to the BBMG jets (for BBMG source feedback)

  //============================================================================
  // Internal helpers
  //============================================================================

  /// @brief Cartesian -> Milne boost of a source 4-position and 4-momentum,
  /// applied in place. No-op unless is_hyperbolic; a no-op too if outside the
  /// light cone (t*t < z*z). Static: it only needs the coordinate flag.
  static void transform_to_hyperbolic( bool is_hyperbolic,
      double& px, double& py, double& pz,
      double& x,  double& y,  double& z,
      double& t,  double& E );

  /// @brief Common time cut for the file readers: skip entries outside the light
  /// cone (hyperbolic) and keep only those whose start time (tau or t) exceeds t0.
  bool passes_time_cut(double t, double z) const;

  /// @brief Allocate + default-initialize the propagating-source master state.
  void allocate_master_state(std::size_t n);

  /// @brief Fill source slot i of the `sources` AoSoA (applying the Milne boost)
  /// and, when seed_master, seed the propagating-source master state from the same
  /// entry. Shared by the AMPT, SMASH and JEWEL readers.
  void register_source(std::size_t i,
                       double px, double py, double pz, double E, double mass,
                       double x, double y, double z, double t,
                       int pid, double baryon_charge, double strangeness_charge,
                       double electric_charge, bool seed_master);

  /// @brief Log one loaded source entry (index i) through formatted_output.
  void log_source_entry(std::size_t i);

  // File readers / real-time builders
  void read_source_file();
  void read_ampt_file();
  void read_smash_file();
  void read_jewel_file();

  void build_jet_sources(double t, double t_prev);
  void build_bbmg_sources(double t, double t_prev);
  void build_toy_sources(double t, double t_prev);
  void build_step_sources(double t, double t_prev);

  /// @brief Lifetime / escape diagnostics for the propagating sources, written to
  /// <results>/smash_source_report.txt.
  void write_source_report();

public:

  /// @brief Constructor for the Source class.
  Source( std::shared_ptr<Settings> settingsPtr_in,
          std::shared_ptr<SystemState<D>> systemPtr_in )
  {
    initialize(settingsPtr_in, systemPtr_in);
  }

  ~Source(){};

  /// @brief Provide the BBMG jet object whose per-step energy loss is deposited
  /// into the fluid when source_model == "BBMG".
  void set_bbmg(BBMG<D>* b) { bbmgPtr = b; }

  void initialize( std::shared_ptr<Settings> settingsPtr_in,
                   std::shared_ptr<SystemState<D>> systemPtr_in );

  /// @brief Deposit all active, in-window source currents into the fluid. Called
  /// once per RK substep from SPHWorkstation::get_time_derivatives.
  void deposit_sources();

  /// @brief Advance + deplete every propagating IC source by ONE timestep. Called
  /// once per timestep from SPHWorkstation::advance_timestep (next to bbmg.propagate()).
  void advance_sources();

};
}// namespace ccake



#endif // SOURCE_H
