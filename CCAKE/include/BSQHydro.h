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

  static constexpr int rk_order = 2; ///TODO: make this a setting
  std::shared_ptr<Settings> settingsPtr; ///< Object containing settings parsed from input file
  std::shared_ptr<SystemState<D>> systemPtr; ///< Object containing the SPH System (linked list, particles, etc.)
  std::shared_ptr<SPHWorkstation<D,TEOM>> wsPtr; ///< Object containing the kernel function and its derivatives
  std::shared_ptr<Output<D>> outPtr; ///< Input/Output object
};
}
#endif
