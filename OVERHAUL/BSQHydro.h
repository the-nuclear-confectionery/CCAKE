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

  void initialize_entropy_and_charge_densities();
  void initial_smoothing();

}

#endif