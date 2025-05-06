#ifndef SPH_WORKSTATION_H
#define SPH_WORKSTATION_H
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>

#include <Cabana_Core.hpp>

#include "constants.h"
//#include "eom.h"
#include "eom_default.h"
#include "eom_cartesian.h"
#include "evolver.h"
#include "formatted_output.h"
#include "freeze_out.h"
#include "kernel.h"
#include "settings.h"
#include "stopwatch.h"
#include "system_state.h"
#include "transport_coefficients.h"
#include "thermodynamic_info.h"
#include "eos_interpolator.h"




namespace ccake
{
/// @class SPHWorkstation	
/// @brief A class to hold the SPH system and perform SPH operations.
/// @details This class holds the SPH system of particles and performs the loops
/// over the particles to calculate the SPH quantities. It also is responsible 
/// for invoking the evolver and handle to it equations of motion. The class is
/// templated on the dimensionality of the simulation and the equation of motion
/// class to be used. To use a different equation of motion, a new class 
/// must be created and should be passed as the `TEOM` template argument.
///
/// The functions the EoM class denoted in this document as (`TEOM`) are:
/// - `TEOM( std::shared_ptr<Settings> settingsPtr_in )`
/// - `KOKKOS_FUNCTION static double gamma_calc(double u[D], const double &time_squared)`
/// - `KOKKOS_FUNCTION static double get_LRF(const double &lab, const double &gamma, const double &time_squared)`
/// - `static void reset_pi_tensor(std::shared_ptr<SystemState<D>> sysPtr)`
/// - `static void evaluate_time_derivatives( std::shared_ptr<SystemState<D>> sysPtr)`
///
/// @note Creating an EoM class that inherits from a base class is not supported.
/// Besides a potential slow down coming from additional memory access that
/// must be done to access virtual functions, in some plataforms (such as GPUs)
/// this is not supported, compromising the portability of the code. For this
/// reason, the template architecture was chosen.
///
/// @tparam D Is the dimensionality of the simulation.
/// @tparam TEOM is the equation of motion class to be used in the simulation.
template <unsigned int D, template<unsigned int> class TEOM>
class SPHWorkstation
{
private:

  static constexpr int    VERBOSE        = 3;
  static constexpr double TOLERANCE      = 0.0;
  static constexpr bool   REGULATE_LOW_T = false;
  double total_energy_loss = 0.0;
  bool hit = false; ///< Flag to check if the source has hit
  std::shared_ptr<Settings> settingsPtr;     ///< Object containing settings parsed from input file
  std::shared_ptr<SystemState<D>> systemPtr; ///< Object containing the SPH System (linked list, particles, etc.)
  // evolver
  Evolver<D> evolver;

  // equation of state
  EquationOfState eos;
  std::shared_ptr<EoS_Interpolator> eos_interpolatorPtr; ///< Object to interpolate the EoS

  transport_coefficients::parameters transp_coeff_params; ///< Transport coefficients parameters


  void add_buffer(double default_e);
public:

  std::shared_ptr<FreezeOut<D>> freezePtr; ///< Object storing the freeze out data
  /// @brief Set up the freeze out object.
  void setup_freeze_out(){
    freezePtr = std::make_shared<FreezeOut<D>>(settingsPtr, systemPtr);}
  void freeze_out_particles();

  //============================================================================
  // Constructors/destructors
  SPHWorkstation() = delete; ///< Default constructor is deleted. System state and settings must be passed in.
  /// @brief Constructor for SPHWorkstation.
  /// @details This constructor initializes the SPHWorkstation object with the
  /// settings and system state objects. It also initializes the evolver object
  /// with the same settings and system state objects.
  /// @param settingsPtr_in A shared pointer to the settings object.
  /// @param systemPtr_in A shared pointer to the system state object.
  SPHWorkstation( shared_ptr<Settings> settingsPtr_in,
                  shared_ptr<SystemState<D>> systemPtr_in )
    : settingsPtr(settingsPtr_in),
      systemPtr(systemPtr_in),
      evolver(Evolver(settingsPtr_in, systemPtr_in)) {}
  KOKKOS_FUNCTION
  ~SPHWorkstation(){}

  //============================================================================
  // initialize workstation (now includes eos initialization)
  void initialize();

  //============================================================================
  // routines for resetting quantities
  void reset_pi_tensor(double time_squared);
  void calculate_gamma_and_velocities();
  //void regulator();
  void process_initial_conditions();
  void initialize_entropy_and_charge_densities();
  void initial_smoothing();

  //============================================================================
  // smoothing
  void calculate_intial_sigma(double time_squared);
  void smooth_all_particle_fields(double time_squared);
  void smooth_all_particle_gradients(double time_squared);

  //============================================================================
  // functions to apply action to all particles
  //-------------------------------------------
  void update_all_particle_thermodynamics();
  void update_all_particle_viscosities();
  void get_time_derivatives();

  //============================================================================
  // routines to edit particles directly
  double locate_phase_diagram_point_eBSQ( Particle<D> & p,
                 double e_In, double rhoB_In, double rhoS_In, double rhoQ_In );
  double locate_phase_diagram_point_eBSQ(Particle<D> & p, double e_In);
  //KOKKOS_INLINE_FUNCTION
  void locate_phase_diagram_point_sBSQ(Particle<D> & p, double s_In,
          double rhoB_In, double rhoS_In, double rhoQ_In );
  //KOKKOS_INLINE_FUNCTION
  void locate_phase_diagram_point_sBSQ(Particle<D> & p, double s_In );

  void set_bulk_Pi();
  void calculate_extensive_shv();
  void set_diffusion();
  bool continue_evolution();
  void advance_timestep( double dt, int rk_order );
  

};

} // namespace ccake
#endif