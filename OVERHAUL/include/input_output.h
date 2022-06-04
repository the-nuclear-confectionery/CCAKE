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
#include "sph_workstation.h"
#include "system_state.h"

// Forward declaration of friend classes
class EquationOfState;
class Settings;
class SPHWorkstation;
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
  void set_SPHWorkstationPtr( SPHWorkstation * wsPtr_in );

  int n_timesteps_output = 0;

  void load_settings_file( string path_to_settings_file );
  void set_EoS_type();
  void set_results_directory( string path_to_results_directory );
  void read_in_initial_conditions();
  void add_buffer();

  void print_system_state();
  void print_system_state_to_txt();
  void print_system_state_to_HDF();
  void print_freeze_out();


  inline void remove_char( std::string & s, char c )
              { s.erase(std::remove(s.begin(), s.end(), c), s.end()); }

  void set_value( setting_map & values, const string key, const string value )
  {
    try { values.at(key) = value; }
    catch (const std::out_of_range& oor)
        { std::cerr << endl << "Invalid key: \"" << key
                    << "\" not found in values!" << endl << endl;
          abort(); }
  }

  string get_value( setting_map & values, const string & key )
  {
    try { return values.at(key); }
    catch (const std::out_of_range& oor)
        { std::cerr << endl << "Invalid key: \"" << key
                    << "\" not found in values!" << endl << endl;
          abort(); }
  }

private:

  const char comment_character = '#'; //make it something else if you want

  string input_directory;
  string output_directory;

  // these allow I/O to access other objects in BSQHydro
  EquationOfState   * eosPtr      = nullptr;
  Settings          * settingsPtr = nullptr;
  SystemState       * systemPtr   = nullptr;
  SPHWorkstation    * wsPtr       = nullptr;

  interface_to_HDF5 hdf5_file;


public:

  void print_conservation_status()
  {
//    out << setw(12) << setprecision(10) << "t="
//        << systemPtr->t      << " " << systemPtr->Eloss  << " "
//        << systemPtr->E0     << " " << systemPtr->Etot   << " "
//        << systemPtr->S      << " " << systemPtr->Btotal << " "
//        << systemPtr->Stotal << " " << systemPtr->Qtotal << endl;

    // print energy/entropy and conserved charge totals
    stringstream ss;
    ss  << "t = "
        << systemPtr->t      << ": " << scientific        << setw(10)
        << systemPtr->Eloss  << " "  << systemPtr->S      << " " 
        << systemPtr->Btotal << " "  << systemPtr->Stotal << " "
        << systemPtr->Qtotal << defaultfloat;
    formatted_output::summarize(ss.str());
  }

};



#endif