#ifndef OUTPUT_H
#define OUTPUT_H

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

#include <Cabana_Core.hpp>

#include "constants.h"
#include "defaults.h"
#include "formatted_output.h"
#include "interface_to_HDF5.h"
#include "mathdef.h"
#include "vector.h"
#include "particle.h"
#include "system_state.h"
#include "vector.h"
#include "freeze_out.h"

namespace ccake{
template<unsigned int D>
class Output
{
public:

  Output(  std::shared_ptr<Settings> settingsPtr_in,
           std::shared_ptr<SystemState<D>> sys_in);
  ~Output();

  void set_results_directory( string path_to_results_directory );
  void print_system_state();
  void print_freeze_out(std::shared_ptr<FreezeOut<D>> freeze_out);
  void print_conservation_status();
  std::string get_freeze_out_filename(){return output_directory + "/freeze_out.dat";}


private:

  int n_timesteps_output = 0;

  string input_directory;
  string output_directory;

  // these allow I/O to access other objects in BSQHydro
  std::shared_ptr<Settings> settingsPtr = nullptr;
  std::shared_ptr<SystemState<D>> systemPtr   = nullptr;

  interface_to_HDF5 hdf5_file;

  void print_system_state_to_txt();
  void print_system_state_to_HDF();

};}


#endif