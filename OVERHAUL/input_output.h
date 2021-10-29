#ifndef INPUT_OUTPUT_H
#define INPUT_OUTPUT_H

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>


#include "mathdef.h"
#include "vector.h"
#include "particle.h"
#include "system_state.h"
#include "equations_of_motion.h"

// Forward declaration of friend classes
class EquationOfState;
class EquationsOfNotion;
class Settings;
class SystemState;

class InputOutput
{
public:

  InputOutput();
  ~InputOutput();

//  void set_EquationOfStatePtr( EquationOfState * eosPtr_in );
  void set_EquationOfStatePtr( std::shared_ptr<EquationOfState> eosPtr_in );
  void set_EquationsOfMotionPtr( EquationsOfMotion * eomPtr_in );
  void set_SettingsPtr( Settings * settingsPtr_in );
  void set_SystemStatePtr( SystemState * systemPtr_in );


  void load_settings_file( string path_to_settings_file ); // load setting
  // paramters for simulation

  void set_EoS_type(); // load in table for interpolation

  void set_results_directory( string path_to_results_directory ); // sets up
  // output directory, will update the outfile as time goes on

  void read_in_initial_conditions(); // talks to
  // system state so that it can set initial system state

  void print_system_state(); //at every time step, will write to output file


private:

  string input_directory;
  string output_directory;

  // these allow I/O to access other objects in BSQHydro
//  EquationOfState   * eosPtr;
    std::shared_ptr<EquationOfState> eosPtr
  EquationsOfMotion * eomPtr;
  Settings          * settingsPtr;
  SystemState       * systemPtr;



};



#endif




















