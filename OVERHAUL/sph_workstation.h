#ifndef SPH_WORKSTATION_H
#define SPH_WORKSTATION_H

#include "system_state.h"
#include "settings.h"
#include "kernel.h"

class SPHWorkstation
{
public:

  SPHWorkstation(){}
  ~SPHWorkstation(){}

  void set_SystemStatePtr( SystemState * systemPtr_in );
  void set_SettingsPtr( Settings * settingsPtr_in );

  void setshear();

  void initialize_entropy_and_charge_densities();
  void initial_smoothing();

  void smooth_fields(int a, bool init_mode = false);
  void smooth_gradients( int a, double tin, int & count );

private:
  
  SystemState * systemPtr;
  Settings * settingsPtr;


}



#endif