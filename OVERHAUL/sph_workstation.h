#ifndef SPH_WORKSTATION_H
#define SPH_WORKSTATION_H

#include "equations_of_motion.h"
#include "kernel.h"
#include "settings.h"
#include "system_state.h"

class SPHWorkstation
{
public:

  SPHWorkstation(){};
  ~SPHWorkstation(){};

  void set_EquationOfStatePtr( EquationOfState * eosPtr_in );
  void set_EquationsOfMotionPtr( EquationsOfMotion * eomPtr_in );
  void set_SystemStatePtr( SystemState * systemPtr_in );
  void set_SettingsPtr( Settings * settingsPtr_in );

  void setshear(bool is_first_timestep);

  void process_initial_conditions();
  void initialize_entropy_and_charge_densities();
  void initial_smoothing();

  void smooth_fields(int a, bool init_mode = false);
  void smooth_gradients( int a, double tin, int & count );

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

private:
  
  SystemState * systemPtr   = nullptr;
  Settings * settingsPtr    = nullptr;
  EquationOfState * eosPtr  = nullptr;
  EquationsOfMotion * eomPtr = nullptr;


};



#endif