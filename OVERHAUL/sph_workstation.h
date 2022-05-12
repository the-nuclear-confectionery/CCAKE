#ifndef SPH_WORKSTATION_H
#define SPH_WORKSTATION_H

#include <omp.h>

#include "eom.h"
#include "eom_default.h"
#include "freeze_out.h"
#include "kernel.h"
#include "settings.h"
#include "Stopwatch.h"
#include "system_state.h"
#include "transport_coefficients.h"

class SPHWorkstation
{
friend class BSQHydro;
friend class InputOutput;

private:
  
  static constexpr int    VERBOSE        = 0;
  static constexpr double TOLERANCE      = 0.0;
  static constexpr bool   REGULATE_LOW_T = false;

  Settings        * settingsPtr          = nullptr;
  SystemState     * systemPtr            = nullptr;

  // equations of motion
  pEquationsOfMotion pEoM;

  // equation of state
  EquationOfState eos;

  // transport coefficients
  TransportCoefficients tc;

  // freeze out
  FreezeOut fo;


public:

  // default constructor/destructor
  SPHWorkstation(){}
  ~SPHWorkstation(){}

  // initialize pointers
  void set_SystemStatePtr( SystemState * systemPtr_in );
  void set_SettingsPtr( Settings * settingsPtr_in );

  // initialize workstation (now includes eos initialization)
  void initialize()
  {
    // set up equation of motion object
    pEoM = std::make_shared<EoM_default>();
    pEoM->set_SettingsPtr( settingsPtr );

    // set up equation of state
    eos.set_SettingsPtr( settingsPtr );
    eos.init();

    // set up transport coefficients
    tc.set_SettingsPtr( settingsPtr );
    tc.initialize();

    // set up freeze out (constant energy density
    fo.set_EquationOfStatePtr( eosPtr );
    fo.set_SettingsPtr( settingsPtr );
    fo.set_SystemStatePtr( systemPtr );
    systemPtr->efcheck = eos.efreeze(settingsPtr->Freeze_Out_Temperature);
    systemPtr->sfcheck = eos.sfreeze(settingsPtr->Freeze_Out_Temperature);
    fo.initialize( systemPtr->efcheck );
  }

  // routines for resetting quantities
  void reset_linklist() { systemPtr->linklist.reset(); }
  void reset_pi_tensor();

  void process_initial_conditions();
  void initialize_entropy_and_charge_densities();
  void initial_smoothing();


  // smoothing
  void smooth_fields( Particle & pa );
  void smooth_gradients( Particle & pa, double tin );

  // functions to apply action to all particles
  void smooth_all_particle_fields()
        { Stopwatch sw;
          sw.Start();
          for ( auto & p : systemPtr->particles )
            smooth_fields(p);
          sw.Stop();
          cout << "t=" << systemPtr->t << ": finished " << __FUNCTION__ << " in " << sw.printTime()
              << " s using " << omp_get_num_threads() << " threads.\n";
        }

  void smooth_all_particle_gradients()
        {
          Stopwatch sw;
          sw.Start();
          #pragma omp parallel for schedule(dynamic)
          for ( auto & p : systemPtr->particles )
            smooth_gradients( p, systemPtr->t );
          sw.Stop();
          cout << "t=" << systemPtr->t << ": finished " << __FUNCTION__ << " in " << sw.printTime()
              << " s using " << omp_get_num_threads() << " threads.\n";
        }

  void update_all_particle_thermodynamics()
        { for ( auto & p : systemPtr->particles )
            calcbsq( p ); }

  void update_all_particle_viscosities()
        { for ( auto & p : systemPtr->particles )
            setvisc( p ); }



  void evaluate_all_particle_time_derivatives()
        { for ( auto & p : systemPtr->particles )
          { p.hydro.ID = p.ID;
            p.hydro.t  = systemPtr->t;
            pEoM->evaluate_time_derivatives( p.hydro, p.thermo );} }


  int do_freezeout_checks();
  void update_freeze_out_lists();
  void finalize_freeze_out(int curfrz);


  // routines to edit particles directly
  double locate_phase_diagram_point_eBSQ(
          Particle & p, double e_In,
          double rhoB_In, double rhoS_In, double rhoQ_In );
  double locate_phase_diagram_point_eBSQ(
          Particle & p, double e_In );
  void locate_phase_diagram_point_sBSQ(
          Particle & p, double s_In,
          double rhoB_In, double rhoS_In, double rhoQ_In );
  void locate_phase_diagram_point_sBSQ(
          Particle & p, double s_In );

  void calcbsq( Particle & p );
  void setvisc( Particle & p );


  // freeze-out routines
//  void bsqsvfreezeout(int curfrz);
//  void bsqsvinterpolate(int curfrz);

  // misc. routine
  double gradPressure_weight(const int a, const int b);



  // Move this into a different namespace or something?
  // It feels like this should be organized separately
  void advance_timestep_rk2( double dt );
  void advance_timestep_rk4( double dt );
  void advance_timestep( double dt, int rk_order )
  {
    // turn on freeze-out flag initially
    systemPtr->cfon = 1;

    switch ( rk_order )
    {
      case 2:
        advance_timestep_rk2( dt );
        break;
      case 4:
        advance_timestep_rk4( dt );
        break;
      default:
        std::cerr << "Invalid Runge-Kutta order!" << std::endl;
        exit(8);
        break;
    }

    // set number of particles which have frozen out
    systemPtr->number_part = systemPtr->get_frozen_out_count();
//    std::cout << "Check termination conditions: "
//              << system.t << "   "
//              << settings.tend << "   "
//              << system.get_frozen_out_count() << "   "
//              << system.number_part << "   "
//              << system.n() << "   "
//              << ( system.t < settings.tend ) << "   "
//              << ( system.number_part < system.n() ) << endl;

    return;
  }

  //MOVE THIS TO ITS OWN CLASS USING TRAVIS' IMPROVEMENTS
  void get_time_derivatives();

  // decide whether to continue evolving
  bool continue_evolution()
  {
    return ( systemPtr->t < settingsPtr->tend )
            && ( systemPtr->number_part < systemPtr->n() );
  }

};



#endif