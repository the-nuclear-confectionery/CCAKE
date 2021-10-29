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

//  void set_EquationOfStatePtr( EquationOfState * eosPtr_in );
  void set_EquationOfStatePtr( std::shared_ptr<EquationOfState> eosPtr_in );
  void set_SystemStatePtr( SystemState * systemPtr_in );
  void set_SettingsPtr( Settings * settingsPtr_in );

  void setshear();

  void process_initial_conditions();
  void initialize_entropy_and_charge_densities();
  void initial_smoothing();

  void smooth_fields(int a, bool init_mode = false);
  void smooth_gradients( int a, double tin, int & count );

private:
  
  SystemState * systemPtr;
  Settings * settingsPtr;
//  EquationOfState * eosPtr;
    std::shared_ptr<EquationOfState> eosPtr = std::make_shared<EquationOfState>();


};



#endif