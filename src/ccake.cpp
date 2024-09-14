//============================================================================//
// Code:    BSQHydro
// Authors: Dekra Almaalol, Travis Dore, Jaki Noronha-Hostler, Lydia Spychalla,
//          Christopher Plumberg, Nikolas Cruz-Camacho
// Date:    October 7, 2021
// Purpose: Run boost-invariant event-by-event hydrodynamics with conserved BSQ
//          charges using the smoothed particle hydrodynamics (SPH) formalism
//============================================================================//
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../include/BSQHydro.h"
#include "../include/formatted_output.h"
#include "../include/welcome.h"

#ifdef DEBUG
#include <fenv.h> 
#endif

using std::cout;
using std::endl;
using std::string;
using std::vector;

int main (int argc, char *argv[])
{
  #ifdef DEBUG
  //feenableexcept(FE_INVALID | FE_OVERFLOW);
  #endif
  // Print the welcome message.
  cout << "Made it to line " << __LINE__ << endl;
  message::welcome();


  // Check if proper command-line arguments passed.
  cout << "Made it to line " << __LINE__ << endl;
  if ( argc < 3 )
  {
    std::cerr << "Usage: ./ccake /path/to/settings/file "
                 "/path/to/results/directory" << std::endl;
    std::cerr << "Please cite all our papers." << std::endl;
    exit(1);
  }

  cout << "Made it to line " << __LINE__ << endl;
  //----------------------------------------------
  formatted_output::announce("Reading in command-line arguments");
  cout << "Made it to line " << __LINE__ << endl;

  // This is where all parameters are initialized.
  cout << "Made it to line " << __LINE__ << endl;
  string path_to_settings_file     = argv[1];
  string path_to_results_directory = argv[2];


  //----------------------------------------------
  cout << "Made it to line " << __LINE__ << endl;
  formatted_output::report( "Input parameters file: "
                            + path_to_settings_file );
  formatted_output::report( "All results will be stored in: "
                            + path_to_results_directory );  


  // Define and set up the simulation object itself.
  cout << "Made it to line " << __LINE__ << endl;
  BSQHydro simulation;
  simulation.set_results_directory( path_to_results_directory );


  //----------------------------------------------
  cout << "Made it to line " << __LINE__ << endl;
  formatted_output::announce("Loading data");


  // Load file containing parameter settings.
  cout << "Made it to line " << __LINE__ << endl;
  simulation.load_settings_file( path_to_settings_file );


  // Read in initial conditions (type/path defined in path_to_settings_file).
  cout << "Made it to line " << __LINE__ << endl;
  simulation.read_in_initial_conditions();


  // This is where the hydrodynamic simulation is set up and initialized.
  cout << "Made it to line " << __LINE__ << endl;
  simulation.initialize_hydrodynamics();


  // Duh.
  cout << "Made it to line " << __LINE__ << endl;
  simulation.run();


  // Print success message.
  formatted_output::announce("Summary: hydrodynamic evolution completed "
                             "successfully");

  return 0;
}
