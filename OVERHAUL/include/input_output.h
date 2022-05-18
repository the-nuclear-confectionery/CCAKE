#ifndef INPUT_OUTPUT_H
#define INPUT_OUTPUT_H

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "defaults.h"
#include "interface_to_HDF5.h"
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
  void print_system_state_to_txt();
  void print_system_state_to_HDF();

  inline void remove_char( std::string & s, char c )
              { s.erase(std::remove(s.begin(), s.end(), c), s.end()); }

  void set_value( setting_map & values, const string key, const string value )
  {
    try { cout << __FUNCTION__ << ": " << key << " " << value << endl;
          auto & tmp = values.at(key);
          tmp = value; }
    catch (const std::out_of_range& oor)
        { std::cerr << endl << "Invalid key: \"" << key
                    << "\" not found in values!" << endl << endl;
          abort(); }
  }

  string get_value( setting_map & values, const string & key )
  {
    try { return values.at(key); }
    catch (const std::out_of_range& oor)
        { std::cerr << "Invalid key: \"" << key << "\" not found in values!\n"; abort(); }
  }

private:

  string input_directory;
  string output_directory;

  // these allow I/O to access other objects in BSQHydro
  EquationOfState   * eosPtr      = nullptr;
  Settings          * settingsPtr = nullptr;
  SystemState       * systemPtr   = nullptr;


  interface_to_HDF5 hdf5_file;


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