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

using std::string;
using std::vector;

namespace ccake{

template<unsigned int D>
class BSQHydro
{

public:

  BSQHydro() = delete;
  BSQHydro( std::shared_ptr<Settings> settingsPtr_in );
  ~BSQHydro(){}

  void read_in_initial_conditions();
  void initialize_hydrodynamics();
  void run();
  void find_freeze_out_surface();
  void print_results();

private:

  static constexpr int rk_order = 4; //TODO: make this a setting

  string input_directory; ///< Path to directory containing input files
  string output_directory; ///< Path to directory where output files will be stored

  std::shared_ptr<Settings> settingsPtr; ///< Object containing settings parsed from input file

  InputOutput<D> io; ///< Input/Output object

  // hold freeze-out surface
  //FreezeOutSurface freeze_out_surface;

  SystemState<D> system; ///< Object containing the SPH System (linked list, particles, etc.)

  // the workstation for performing SPH-related actions on the system
  SPHWorkstation<D> ws;

};
}
#endif
