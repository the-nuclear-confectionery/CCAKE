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
#include "../include/welcome.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

int main (int argc, char *argv[])
{
  // Print the welcome message.
  message::welcome();
  //message::guilty_welcome();

  // Check if proper command-line arguments passed.
  if ( argc < 3 )
  {
    std::cerr << "Usage: ./persephone /path/to/settings/file "
                 "/path/to/results/directory" << std::endl;
    std::cerr << "Please cite all our papers." << std::endl;
    exit(1);
  }

  // This is where all parameters are initialized.
  string path_to_settings_file     = argv[1];
  string path_to_results_directory = argv[2];

  // Define and set up the simulation object itself.
  BSQHydro simulation;
  simulation.set_results_directory( path_to_results_directory );
  simulation.load_settings_file( path_to_settings_file );

  // Read in initial conditions (type/path defined in path_to_settings_file).
  simulation.read_in_initial_conditions();

  // This is where the hydrodynamic simulation is set up and initialized.
  simulation.initialize_hydrodynamics();

  // Duh.
  simulation.run();

  // Construct freeze-out surface.
  simulation.find_freeze_out_surface();

  // Save any relevant output.
  simulation.print_results();

  // Farewell message.
  std::cout << "Simulation ran successfully!" << std::endl;

  return 0;
}
