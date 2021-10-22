#ifndef IO_H
#define IO_H

#include "constants.h"

using std::string;

class InputOutput
{

public:

  // input functions
  void read_in_settings_file( string settings_filename );
  void read_in_initial_conditions_file( string initial_conditions_filename );

  // output functions
  void print_simulation_settings();
  void print_hydrodynamic_state( const System_state & system );
  void print_freezeout_surface( const Freeze_out_surface & surface );

private:

}

#endif