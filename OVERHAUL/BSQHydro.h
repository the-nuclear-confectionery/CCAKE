#ifndef BSQHYDRO_H
#define BSQHYDRO_H

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

#include "constants.h"
#include "eos.h"
#include "equations_of_motion.h"
#include "input_output.h"
#include "system_state.h"
#include "settings.h"
#include "kernel.h"
#include "sph_workstation.h"

using std::string;
using std::vector;

class BSQHydro
{

public:

  BSQHydro(){}
  ~BSQHydro(){}

  void load_settings_file( string path_to_settings_file );
  void set_results_directory( string path_to_results_directory );
  void read_in_initial_conditions();
  void initialize_hydrodynamics();
  void run();
  void find_freeze_out_surface();
  void print_results();



private:

  string input_directory;
  string output_directory;


  // equation of state information
  EquationOfState eos;

  // equation of motion object
  EquationsOfMotion eom;

  // input/output
  InputOutput io;

  // hold freeze-out surface
  //FreezeOutSurface freeze_out_surface;

  // the current state of the hydrodynamic simulation
  SystemState system;

  // the workstation for performing SPH-related actions on the system
  SPHWorkstation ws;

  // all settings for the hydro simulation
  Settings settings;


  struct Input_Parameters
  {
    string IC_type; // specify initial condition type
    double h; // static SPH cutoff paramter
    double dt; // time step in fm
    double t0; // initial time in fm
    string EoS_type; // specify equation of state type
    string EoS_option; // specify specifc option for EOS
    // there should an associated EoS directory with tables
    string eta; // specificy the shear viscosity type to use
    // in transport cpefficient file
    string zeta; // specificy the bulk viscosity type to use
    // in transport cpefficient file
    double Freeze_Out_Temperature;
    string Freeze_Out_Type;
  };

  struct Initial_Conditions
  {
    vector<string> headers; // vector of header parameters as string
    vector<vector<double> > density_grid; // the hydro grid, read in and stored
  };

  Input_Parameters input_parameters;
  Initial_Conditions initial_conditions;

  void trim_initial_conditions();
};

#endif
