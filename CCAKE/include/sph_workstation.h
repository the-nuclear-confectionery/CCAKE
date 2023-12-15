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
#include "evolver.h"
#include "formatted_output.h"
#include "freeze_out.h"
#include "kernel.h"
#include "settings.h"
#include "stopwatch.h"
#include "system_state.h"
#include "transport_coefficients.h"
#include "thermodynamic_info.h"
//#define ONLINE_INVERTER
#ifndef ONLINE_INVERTER
#include "eos_interpolator.h"
#endif




namespace ccake
{

/// @brief A class to hold the SPH system and perform SPH operations.
/// @details This class holds the SPH system of particles and performs the loops
/// over the particles to calculate the SPH quantities. It also contains the
/// equations of motion and the evolver. The class is templated on the
/// dimensionality of the simulation and the equation of motion class to be
/// used. To use a different equation of motion, the user must create a new
/// class similar to the default equation of motion class (eom_default.h) and
/// pass it as the second template argument to the SPHWorkstation class.
///
/// Creating an EoM class that inherits from the default is not supported and
/// disencouraged, since it will require in each function call, at runtime,
/// for the process to decide which equation of motion to use, potentially
/// slowing down the code. It's the user responsibility to make sure that their
/// own EoM class implements all the functions required by the SPHWorkstation.
///
/// The functions the Eom class must implement are:
/// - constructor: must receive a pointer to the Settings object as input
///
/// @tparam TEOM is the equation of motion class to be used in the simulation.
/// @tparam D Is the dimensionality of the simulation.
template <unsigned int D, template<unsigned int> class TEOM>
class SPHWorkstation
{
private:

  static constexpr int    VERBOSE        = 3;
  static constexpr double TOLERANCE      = 0.0;
  static constexpr bool   REGULATE_LOW_T = false;

  std::shared_ptr<Settings> settingsPtr; ///< Object containing settings parsed from input file
  std::shared_ptr<SystemState<D>> systemPtr; ///< Object containing the SPH System (linked list, particles, etc.)
  // evolver
  Evolver<D> evolver;

  // equation of state
  EquationOfState eos;
  #ifndef ONLINE_INVERTER
  std::shared_ptr<EoS_Interpolator> eos_interpolatorPtr;
  #endif

  // transport coefficients parameters
  transport_coefficients::parameters transp_coeff_params;


  void add_buffer(double default_e);
public:

  // freeze out
  std::shared_ptr<FreezeOut<D>> freezePtr;
  void setup_freeze_out(){freezePtr = std::make_shared<FreezeOut<D>>(settingsPtr, systemPtr);}
  void freeze_out_particles();

  //============================================================================
  // Constructors/destructors
  SPHWorkstation() = delete; ///< Default constructor is deleted. System state and settings must be passed in.
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

  void process_initial_conditions();
  void initialize_entropy_and_charge_densities();
  void initial_smoothing();

  //============================================================================
  // smoothing
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
  bool continue_evolution();
  void advance_timestep( double dt, int rk_order );

};

} // namespace ccake
#endif