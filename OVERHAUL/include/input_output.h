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

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

// Forward declaration of friend classes
class EquationOfState;
class Settings;
class SystemState;

class InputOutput
{
public:

  InputOutput();
  ~InputOutput();

  int initialize_HDF();

  void set_EquationOfStatePtr( EquationOfState * eosPtr_in );
  void set_SettingsPtr( Settings * settingsPtr_in );
  void set_SystemStatePtr( SystemState * systemPtr_in );

  int n_timesteps_output = 0;

  void load_settings_file( string path_to_settings_file );
  void set_EoS_type();
  void set_results_directory( string path_to_results_directory );
  void read_in_initial_conditions();

  void print_system_state();
  void print_shear();

private:

  string input_directory;
  string output_directory;

  // these allow I/O to access other objects in BSQHydro
  EquationOfState   * eosPtr      = nullptr;
  Settings          * settingsPtr = nullptr;
  SystemState       * systemPtr   = nullptr;

  H5File file;
 	string GROUPEVENT_NAME = "/Event";

  string get_zero_padded_int( int i, int width );
  void output_double_attribute( Group & group, double value, string name );
  void output_dataset( string FRAME_NAME, const double time );
  void set_units(DataSet & ds, const std::string & units);



public:

  void print_conservation_status(std::ostream & out = std::cout)
  {
    // print energy/entropy and conserved charge totals
    out << setw(12) << setprecision(10) << "t="
        << systemPtr->t      << " " << systemPtr->Eloss  << " "
        << systemPtr->E0     << " " << systemPtr->Etot   << " "
        << systemPtr->S      << " " << systemPtr->Btotal << " "
        << systemPtr->Stotal << " " << systemPtr->Qtotal << endl;
  }

};



#endif