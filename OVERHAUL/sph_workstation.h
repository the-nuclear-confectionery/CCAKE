#ifndef SPH_WORKSTATION_H
#define SPH_WORKSTATION_H

#include "kernel.h"
#include "settings.h"
#include "system_state.h"

class SPHWorkstation
{
private:
  
  static constexpr int    VERBOSE        = 0;
  static constexpr double TOLERANCE      = 0.0;
  static constexpr bool   REGULATE_LOW_T = false;

  SystemState     * systemPtr            = nullptr;
  Settings        * settingsPtr          = nullptr;
  EquationOfState * eosPtr               = nullptr;

public:

  // default constructor/destructor
  SPHWorkstation(){}
  ~SPHWorkstation(){}

  // initialize pointers
  void set_EquationOfStatePtr( EquationOfState * eosPtr_in );
  void set_SystemStatePtr( SystemState * systemPtr_in );
  void set_SettingsPtr( Settings * settingsPtr_in );

  // routines for resetting quantities
  void reset_linklist() { systemPtr->linklist.reset(); }
  void reset_pi_tensor();

  void process_initial_conditions();
  void initialize_entropy_and_charge_densities();
  void initial_smoothing();

  void smooth_fields( Particle & pa );
  void smooth_gradients( Particle & pa, double tin, int & count );
  void smooth_all_particle_fields()
        { for ( auto & p : systemPtr->particles )
            smooth_fields(p); }
  void smooth_all_particle_gradients( int & count )
        { for ( auto & p : systemPtr->particles )
            smooth_gradients( p, systemPtr->t, count ); }

  void update_all_particle_thermodynamics()
        { for ( auto & p : systemPtr->particles )
            p.calcbsq( systemPtr->t ); }

  void update_all_particle_viscosities()
        { for ( auto & p : systemPtr->particles )
            p.setvisc( systemPtr->etaconst, systemPtr->bvf, systemPtr->svf,
                       systemPtr->zTc,      systemPtr->sTc, systemPtr->zwidth,
                       systemPtr->visc ); }


  int do_freeze_out_checks();
  void update_all_particle_dsigma_dt();
  void update_all_particle_fluid_quantities()
        { for ( auto & p : systemPtr->particles )
            p.update_fluid_quantities( systemPtr->t ); }

  void update_freeze_out_lists();




  // Move this into a different namespace or something?
  // It feels like this should be organized separately
  void advance_timestep_rk2( double dt );
  void advance_timestep_rk4( double dt );
  void advance_timestep( double dt, int rk_order )
  {
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
    return;
  }

  //MOVE THIS TO ITS OWN CLASS USING TRAVIS' IMPROVEMENTS
  void get_time_derivatives();

};



#endif