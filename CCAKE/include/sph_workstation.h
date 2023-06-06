#ifndef SPH_WORKSTATION_H
#define SPH_WORKSTATION_H
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>

#include "constants.h"
//#include "eom.h"
#include "eom_default.h"
#include "evolver.h"
#include "formatted_output.h"
#include "freeze_out.h"
//#include "kernel.h"
#include "settings.h"
#include "stopwatch.h"
#include "system_state.h"
#include "transport_coefficients.h"



//forward declare template class
//namespace ccake{
//template <unsigned int D>
//class BSQHydro;
//}
//template <unsigned int D> class Evolver;
//template <unsigned int D> class InputOutput;

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
  std::shared_ptr<TEOM<D>> EoMPtr; ///< Equations of motion object

  // evolver
  Evolver<D> evolver;

  // equation of state
  EquationOfState eos;

  // transport coefficients
  TransportCoefficients transport_coefficients;

  // freeze out
  FreezeOut<D> freeze_out;

  void add_buffer(double default_e);

public:

  //============================================================================
  // Constructors/destructors
  SPHWorkstation() = delete; ///< Default constructor is deleted. System state and settings must be passed in.
  SPHWorkstation( shared_ptr<Settings> settingsPtr_in,
                  shared_ptr<SystemState<D>> systemPtr_in )
    : settingsPtr(settingsPtr_in), systemPtr(systemPtr_in) {}
  ~SPHWorkstation(){}

  //============================================================================
  // initialize pointers
  //void set_SystemStatePtr( SystemState<D> * systemPtr_in );
  //void set_SettingsPtr( Settings * settingsPtr_in );

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
  //void smooth_gradients( Particle<D> & pa, double tin );

  //void get_time_derivatives();


  //============================================================================
  // functions to apply action to all particles
  //-------------------------------------------
  
  //-------------------------------------------
  //void smooth_all_particle_gradients()
  //      {
  //        Stopwatch sw;
  //        sw.Start();
  //        for ( auto & p : systemPtr->particles )
  //          smooth_gradients( p, systemPtr->t );
  //        sw.Stop();
  //        formatted_output::update("finished smoothing particle gradients in "
  //                                  + to_string(sw.printTime()) + " s.");
  //      }

  //-------------------------------------------
  //void update_all_particle_thermodynamics()
  //      {
  //        Stopwatch sw;
  //        sw.Start();
  //        for ( auto & p : systemPtr->particles )
  //          set_phase_diagram_point( p );
  //        sw.Stop();
  //        formatted_output::update("got particle thermodynamics in "
  //                                  + to_string(sw.printTime()) + " s."); }

  //-------------------------------------------
  //void update_all_particle_viscosities()
  //      { for ( auto & p : systemPtr->particles )
  //          set_transport_coefficients( p ); }

  //-------------------------------------------
  //void evaluate_all_particle_time_derivatives()
  //      {
  //        Stopwatch sw;
  //        sw.Start();
  //        for ( auto & p : systemPtr->particles )
  //        {
  //          p.hydro.ID = p.ID;
  //          p.hydro.t  = systemPtr->t;
  //          EoMPtr->evaluate_time_derivatives( p.hydro, p.thermo, p.d_dt_spec );
  //        }
  //        sw.Stop();
  //        formatted_output::update("set particle time derivatives in "
  //                                  + to_string(sw.printTime()) + " s."); }


  //-------------------------------------------
  //void freeze_out_particles();


  //============================================================================
  // routines to edit particles directly
  double locate_phase_diagram_point_eBSQ(double e_In,
          double rhoB_In, double rhoS_In, double rhoQ_In );
  double locate_phase_diagram_point_eBSQ(double e_In );
  void locate_phase_diagram_point_sBSQ(Particle<D> & p, double s_In,
          double rhoB_In, double rhoS_In, double rhoQ_In );
  void locate_phase_diagram_point_sBSQ(
          Particle<D> & p, double s_In );

  //void set_phase_diagram_point( Particle & p );
  //void set_transport_coefficients( Particle & p );
//
  void set_bulk_Pi();
//
  //// misc. routine
  //double gradPressure_weight(const int a, const int b);


  //============================================================================
  // what it says on the label
//  void advance_timestep( double dt, int rk_order )
//  {
//    Stopwatch sw;
//    sw.Start();
//
//    // turn on freeze-out flag initially
//    systemPtr->do_freeze_out = true;
//
//    // use evolver to actually do RK evolution
//    // (pass workstation's own time derivatives function as lambda)
//    evolver.advance_timestep( dt, rk_order,
//                              [this]{ this->get_time_derivatives(); } );
//
//    // set number of particles which have frozen out
//    systemPtr->number_part = systemPtr->get_frozen_out_count();
//
//    // keep track of how many timesteps have elapsed
//    systemPtr->number_of_elapsed_timesteps++;
//
//    sw.Stop();
//    formatted_output::report("finished timestep in "
//                              + to_string(sw.printTime()) + " s");
//
//    if ( settingsPtr->max_number_of_timesteps >= 0
//          && systemPtr->number_of_elapsed_timesteps
//              > settingsPtr->max_number_of_timesteps )
//      exit(-1);
//
//    return;
//  }

  //============================================================================
  // decide whether to continue evolving
//  bool continue_evolution()
//  {
//    return ( systemPtr->t < settingsPtr->tend )
//            && ( systemPtr->number_part < systemPtr->n() );
//  }
};

} // namespace ccake


#endif