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

#include "input.h"
#include "BSQHydro.h"
#include "formatted_output.h"
#include "welcome.h"

#ifdef DEBUG
#include <fenv.h>
#endif

namespace cc = ccake;

template<unsigned int D>
void execute_tasks(std::shared_ptr<Settings> settingsPtr)
{
  // Define and set up the simulation object itself.
  cc::BSQHydro<D,cc::EoM_default> simulation(settingsPtr); //TODO: make EoM a setting
                                                       // and a switch case selection.

  formatted_output::announce("Loading initial condition.");
  // Read in initial conditions (type/path defined in path_to_settings_file).
  simulation.read_in_initial_conditions();
  // This is where the hydrodynamic simulation is set up and initialized.;
  simulation.initialize_hydrodynamics();
  //formatted_output::announce("Running hydrodynamics.");
  //Executes the main loop of the simulation.
  simulation.run();

}

int main (int argc, char *argv[])
{
  #ifdef DEBUG
  //feenableexcept(FE_INVALID | FE_OVERFLOW);
  #endif
  // Print the welcome message.
  Kokkos::initialize(argc, argv);
  message::welcome();
  //----------------------------------------------
  formatted_output::announce("Reading in command-line arguments");
  // Read input arguments and parse config file.
  ccake::Input in(argc, argv);

  unsigned int D = in.settingsPtr->dim;
  switch (D)
  {
    case 1:
      formatted_output::announce("Running 1D hydrodynamics");
      execute_tasks<1>(in.settingsPtr);
      break;
    case 2:
      formatted_output::announce("Running 2D hydrodynamics");
      execute_tasks<2>(in.settingsPtr);
      break;
    case 3:
      formatted_output::announce("Running 3D hydrodynamics");
      execute_tasks<3>(in.settingsPtr);
      break;
    default:
      formatted_output::announce("Invalid number of dimensions");
      exit(EXIT_FAILURE);
  }
  // Print success message.
  formatted_output::announce("Summary: hydrodynamic evolution completed "
                             "successfully");
  Kokkos::finalize();

  return 0;
}
