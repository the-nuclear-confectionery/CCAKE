#ifndef BSQHYDRO_H
#define BSQHYDRO_H

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

#include "constants.h"
#include "output.h"
#include "system_state.h"
#include "settings.h"
#include "kernel.h"
#include "sph_workstation.h"
#include "formatted_output.h"
#include "stopwatch.h"
#include "eom_default.h"

using std::string;
using std::vector;

namespace ccake{
/// @class BSQHydro
/// @brief Shell class responsible for setting up and running the hydrodynamic
/// simulation.
///
/// This class is responsible for setting up and running the hydrodynamic
/// simulation. It is the main class that is called from the main function of
/// the code. It is responsible for reading in the initial conditions, setting
/// up the hydrodynamic simulation, and running the main loop of the simulation.
///
/// @todo Reader functions, such as BSQHydro::read_ICCING and 
/// BSQHydro::read_ccake, should be moved to the Input class.
/// @tparam D Dimensionality of the simulation
/// @tparam TEOM Template for the equations of motion to be used
template<unsigned int D,template<unsigned int> class TEOM>
class BSQHydro
{
public:

  BSQHydro() = delete;
  BSQHydro( std::shared_ptr<Settings> settingsPtr_in );
  ~BSQHydro(){}

  void read_in_initial_conditions();
  void initialize_hydrodynamics();
  void run();

private:

  //Auxiliary functions to read ICs
  void read_ICCING();
  void read_ccake();
  void read_dynamical();

  std::shared_ptr<Settings> settingsPtr; ///< Settings object containing configuration parsed from input file
  std::shared_ptr<SystemState<D>> systemPtr; ///< SystemState object containing the SPH System (linked list, particles, etc.)
  std::shared_ptr<SPHWorkstation<D,TEOM>> wsPtr; ///< SPHWorkstation object with the functions executed in the main loop
  std::shared_ptr<Output<D>> outPtr; ///< Output object
};
}
#endif
