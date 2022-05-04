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

#include "mathdef.h"
#include "vector.h"
#include "particle.h"
#include "system_state.h"

// Forward declaration of friend classes
class EquationOfState;
class Settings;
class SystemState;

class InputOutput
{
public:

  InputOutput();
  ~InputOutput();

  void set_EquationOfStatePtr( EquationOfState * eosPtr_in );
  void set_SettingsPtr( Settings * settingsPtr_in );
  void set_SystemStatePtr( SystemState * systemPtr_in );

  int n_timesteps_output = 0;

  void load_settings_file( string path_to_settings_file ); // load setting
  // paramters for simulation

  void set_EoS_type(); // load in table for interpolation

  void set_results_directory( string path_to_results_directory ); // sets up
  // output directory, will update the outfile as time goes on

  void read_in_initial_conditions(); // talks to
  // system state so that it can set initial system state

  void print_system_state(); //at every time step, will write to output file
  void print_shear();        //at every time step, will write to output file

private:

  string input_directory;
  string output_directory;

  // these allow I/O to access other objects in BSQHydro
  EquationOfState   * eosPtr      = nullptr;
  Settings          * settingsPtr = nullptr;
  SystemState       * systemPtr   = nullptr;



public:

  void print_conservation_status(std::ostream & out = std::cout)
  {
    // print energy/entropy and conserved charge totals
    out << setw(12) << setprecision(10)
        << "t=" << systemPtr->t << " "
        << systemPtr->Eloss << " " << systemPtr->E0 << " " << Etot << " "
        << S << " " << Btotal << " " << Stotal << " " << Qtotal << endl;
  }

};



#endif




















