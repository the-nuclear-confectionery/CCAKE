#ifndef SOURCE_H
#define SOURCE_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
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
  /// One entry of the per-step deposit list (this step's energy/momentum chunk).
  using SRC = Cabana::MemberTypes<double[3], // position
                                  double[3], // momentum
                                  double,    // energy
                                  double,    // baryon charge
                                  double,    // strangeness charge
                                  double,    // electric charge
                                  double>;   // renormalization: kernel divisor (1 = none)

  namespace SRC_enum{
    enum SRC_members
    {
      position = 0,
      momentum,
      energy,
      baryon_charge,
      strangeness_charge,
      electric_charge,
      renormalization
    };
  }

  using source_array = Cabana::AoSoA<SRC, DeviceType,VECTOR_LENGTH>;

  #define SRC_VIEW(prefix, SRC_aosoa) \
  auto CONCAT(prefix, position) = Cabana::slice<SRC_enum::position>(SRC_aosoa); \
  auto CONCAT(prefix, momentum) = Cabana::slice<SRC_enum::momentum>(SRC_aosoa); \
  auto CONCAT(prefix, energy) = Cabana::slice<SRC_enum::energy>(SRC_aosoa); \
  auto CONCAT(prefix, baryon_charge) = Cabana::slice<SRC_enum::baryon_charge>(SRC_aosoa); \
  auto CONCAT(prefix, strangeness_charge) = Cabana::slice<SRC_enum::strangeness_charge>(SRC_aosoa); \
  auto CONCAT(prefix, electric_charge) = Cabana::slice<SRC_enum::electric_charge>(SRC_aosoa); \
  auto CONCAT(prefix, renormalization) = Cabana::slice<SRC_enum::renormalization>(SRC_aosoa); \



/// @brief Source terms deposited into the fluid: IC sources read from file
/// (AMPT/SMASH/3dMCGlauber; each releases an energy budget after its formation
/// time) or the per-step energy loss of the BBMG jets.
template <unsigned int D>
class Source
{
private:

  /// One source as read from file (GeV, Cartesian fm), before the Milne transform.
  struct SourceRecord
  {
    double t, x, y, z;
    double E, px, py, pz, mass;
    int    pid;
    double baryon_charge, strangeness_charge, electric_charge;
    int    type = -1, parent = -1;   ///< 3dMCGlauber tags (-1 = untagged)
  };

  /// Deposit on a frozen-at-birth (Freeze = 5) particle; lost only if it never reheats.
  struct PendingFrozenDeposit
  {
    double energy       = 0.0;
    double charge_share = 0.0;
  };

  /// State of one IC source. Units: fm^-1 for energies and momenta.
  struct ICSource
  {
    // read once
    double px, py, pz;        ///< conserved Cartesian momentum
    double mass;
    double E0;                ///< launch energy (Milne p^tau) = initial budget
    double baryon, strange, charge; ///< charges released in one shot on absorption
    double formation_time;    ///< tau (hyperbolic) or t at which the source activates
    int    PID;
    int    type = -1, parent = -1;   ///< 3dMCGlauber element type / parent string (-1 = untagged)

    // evolving
    double x, y, z;           ///< current position (x, y, eta_s in hyperbolic mode)
    double remaining;         ///< energy budget left
    bool   charge_dumped = false;

    // this step (set once per step, read by every RK substep)
    double dE = 0.0;          ///< energy chunk deposited this step
    double vx = 0.0, vy = 0.0, veta = 0.0; ///< Milne velocity (momentum chunk = dE * v)
    bool   dump_now = false;  ///< release the charge this step (absorption)

    // diagnostics
    int    activate_step = -1, absorb_step = -1;  ///< step indices (-1 = never)
    int    defer_steps   = 0;  ///< steps spent waiting for the fluid (kernel_renormalize)
    // this step's kernel coverage (kernel_renormalize only, see compute_coverage)
    double coverage = 1.0;     ///< Shepard coverage sum (~1 where the fluid covers the kernel)
    double frozen_weight = 0.0, permanently_frozen_weight = 0.0;   ///< same sum over Freeze >= 4, = 4
    std::vector<int>    frozen_at_birth_particles;  ///< Freeze = 5 particles under the kernel
    std::vector<double> frozen_at_birth_weights;    ///< ... and their weights

    // deposits onto frozen matter (see record_frozen_deposit)
    double energy_onto_frozen = 0.0;                   ///< delivered energy that landed on Freeze >= 4
    double frozen_share_at_absorb = -1.0;              ///< frozen share when the charge was released
    double energy_onto_permanently_frozen = 0.0;       ///< energy landed on Freeze = 4 (lost)
    double charge_share_onto_permanently_frozen = 0.0; ///< charge share landed on Freeze = 4 (lost)
    std::unordered_map<int, PendingFrozenDeposit> pending_frozen_deposits; ///< particle index -> deposit
  };

  bool has_bbmg_sources = false;
  bool has_ic_sources   = false; ///< AMPT/SMASH/MCGLB sources read from file

  std::vector<ICSource> ic_sources;
  int src_step = 0;   ///< advance_sources call counter

  std::size_t n_sources = 0; ///< Number of entries in the deposit list

  source_array sources; ///< Per-step deposit list

  std::shared_ptr<Settings> settingsPtr    = nullptr; ///< Settings parsed from input file
  std::shared_ptr<SystemState<D>> systemPtr = nullptr; ///< SPH system (linked list, particles, ...)

  BBMG<D>* bbmgPtr = nullptr; ///< Non-owning ref to the BBMG jets

  //============================================================================
  // Internal helpers
  //============================================================================

  /// @brief Cartesian -> Milne transform of a 4-position and 4-momentum, in place
  /// (no-op unless is_hyperbolic, or outside the light cone).
  static void transform_to_hyperbolic( bool is_hyperbolic,
      double& px, double& py, double& pz,
      double& x,  double& y,  double& z,
      double& t,  double& E );

  /// @brief Keep entries inside the light cone whose start time (tau or t) exceeds t0.
  bool passes_time_cut(double t, double z) const;

  /// @brief Turn the records into IC sources.
  void load_sources(const std::vector<SourceRecord>& records);

  /// @brief Fraction of a presented chunk that lands with the capped Shepard factor.
  double delivered_fraction(double coverage) const;

  // File readers / per-step builders
  void read_ampt_file(std::vector<SourceRecord>& records) const;
  void read_smash_file(std::vector<SourceRecord>& records) const;
  void read_mcglb_file(std::vector<SourceRecord>& records) const;

  void build_bbmg_sources(double t);
  void build_step_sources();
  void compute_coverage();
  /// @brief Book the part of this step's deposit that lands on frozen particles.
  void record_frozen_deposit(ICSource& source);

  /// @brief Write <results>/source_ledger.dat (the only source output).
  void write_source_ledger();

public:

  /// @brief Constructor for the Source class.
  Source( std::shared_ptr<Settings> settingsPtr_in,
          std::shared_ptr<SystemState<D>> systemPtr_in )
  {
    initialize(settingsPtr_in, systemPtr_in);
  }

  ~Source(){};

  /// @brief BBMG jets whose energy loss is deposited (source_model == "BBMG").
  void set_bbmg(BBMG<D>* b) { bbmgPtr = b; }

  void initialize( std::shared_ptr<Settings> settingsPtr_in,
                   std::shared_ptr<SystemState<D>> systemPtr_in );

  /// @brief Deposit this step's source currents (once per RK substep).
  void deposit_sources();

  /// @brief Advance and deplete the IC sources (once per timestep).
  void advance_sources();

};
}// namespace ccake



#endif // SOURCE_H
