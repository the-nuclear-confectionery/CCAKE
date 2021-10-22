#include "BSQHydro.h"


// Constructors and destructors.
  BSQHydro::BSQHydro(){}
 ~BSQHydro::BSQHydro(){}


void BSQHydro::load_settings_file( string path_to_settings_file )
{
  io.load_settings_file(path_to_settings_file) //sets the settings path in input_output
  return
}



void BSQHydro::set_results_directory( string path_to_results_directory )
{
  io.set_results_directory(path_to_results_directory) //set the results directory in input_output
  return
}



void BSQHydro::read_in_initial_conditions()
{
  io.read_in_initial_conditions() // tells input_output to talk to system state and set initial system state
  return
}



void BSQHydro::initialize_hydrodynamics()
{
  //chris will write this
}


void BSQHydro::run(){}


void BSQHydro::find_freeze_out_surface(){}


void BSQHydro::print_results(){}