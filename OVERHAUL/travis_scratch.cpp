///////////INPUT_OUTPUT.CPP BEGIN/////////////

#ifndef INPUT_OUTPUT_CPP

load_settings_file( string path_to_settings_file )
{

}


///////////INPUT_OUTPUT.CPP END//////////////


///////////BSQHYDRO.CPP BEGIN/////////////

#include "BSQHydro.h"
#include "input_output.h"


// Constructors and destructors.
  BSQHydro::BSQHydro(){}
 ~BSQHydro::BSQHydro(){}


void BSQHydro::load_settings_file( string path_to_settings_file ){}



void BSQHydro::set_results_directory( string path_to_results_directory ){}



void BSQHydro::read_in_initial_conditions(){}



void BSQHydro::initialize_hydrodynamics(){}


void BSQHydro::run(){}


void BSQHydro::find_freeze_out_surface(){}


void BSQHydro::print_results(){}
///////////BSQHYDRO.CPP END//////////////





///////////BSQHYDRO.H BEGIN/////////////

#ifndef BSQHYDRO_H
#define BSQHYDRO_H

using std::string;

class BSQHydro
{

public:

  BSQHydro();
  ~BSQHydro();

  void load_settings_file( string path_to_settings_file );
  void set_results_directory( string path_to_results_directory );
  void read_in_initial_conditions();
  void initialize_hydrodynamics();
  void run();
  void find_freeze_out_surface();
  void print_results();

  struct Input_Parameters
  {
    string IC_type; // specify initial condition type
    double h; // static SPH cutoff paramter
    double dt; // time step in fm
    double t0; // initial time in fm

  }

  struct Initial_Conditions
  {
    
  }



private:

  string input_directory
  string output_directory;

  // equation of state information
  EquationOfState eos;

  // input/output
  InputOutput io;

  // hold freeze-out surface
  FreezeOutSurface freeze_out_surface;

  // the current state of the hydrodynamic simulation
  SystemState system;

}

#endif

///////////BSQHYDRO.H END//////////////
