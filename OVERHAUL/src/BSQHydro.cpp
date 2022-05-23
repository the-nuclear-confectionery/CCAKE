#include "../include/BSQHydro.h"
#include "../include/formatted_output.h"
#include "../include/settings.h"
#include "../include/input_output.h"
#include "../include/stopwatch.h"

// Constructors and destructors.
BSQHydro::BSQHydro()
{

  // initialize I/O pointers
  io.set_EquationOfStatePtr( &ws.eos );
  io.set_SettingsPtr( &settings );
  io.set_SPHWorkstationPtr( &ws );  // this is probably unnecessary
  io.set_SystemStatePtr( &system );

  // initialize SPH workstation
  ws.set_SystemStatePtr( &system );
  ws.set_SettingsPtr( &settings );
  
  // initialize system state
  system.set_SettingsPtr( &settings );
  
  return;
}

void BSQHydro::load_settings_file( string path_to_settings_file )
{
  // sets the settings path in InputOutput,
  // then loads parameters into Input_parameters struct
  io.load_settings_file(path_to_settings_file);

  // InputOutput talks to EoS and tells it where to find its tables
  io.set_EoS_type();

  return;
}



void BSQHydro::set_results_directory( string path_to_results_directory )
{
  // set the results directory in InputOutput
  io.set_results_directory(path_to_results_directory);
  return;
}



void BSQHydro::read_in_initial_conditions()
{
  // tells InputOutput to talk to system state and set initial system state
  io.read_in_initial_conditions();
  return;
}




////////////////////////////////////////////////////////////////////////////////
void BSQHydro::initialize_hydrodynamics()
{
  formatted_output::announce("Initializing hydrodynamics");
  Stopwatch sw;
  sw.Start();

  // initialize equation of state
  ws.initialize();

  // initialize system state
  system.initialize();

  // trim initial conditions with low-energy density cut-off,
  // filling out initial conditions, and imposing initial freeze-out
  ws.process_initial_conditions();

  // comes after energy cut-off imposed
  system.initialize_linklist();

  // for each particle, find location in phase diagram
  ws.initialize_entropy_and_charge_densities();

  // implement initial smoothing required by SPH formalism
  ws.initial_smoothing();

  // if initializing from full Tmunu, absorb non-equilibrium
  // pressure correction into bulk viscous pressure Pi
  ws.set_bulk_Pi();

  sw.Stop();
  formatted_output::report("hydrodynamics initialization finished in "
                              + to_string(sw.printTime()) + " s");

  return;
}


////////////////////////////////////////////////////////////////////////////////
void BSQHydro::run()
{
  formatted_output::announce("Beginning hydrodynamic evolution");
  Stopwatch sw;
  sw.Start();

  //===================================
  // initialize conserved quantities
  system.conservation_entropy();
  system.conservation_BSQ();

  //===================================
  // print initialized system and status
  io.print_conservation_status();
  io.print_system_state();

  //===================================
  // evolve until simulation terminates
  while ( ws.continue_evolution() )
  {
    //===================================
    // workstation advances by given
    // timestep at given RK order
    ws.advance_timestep( settings.dt, rk_order );

    //===================================
    // re-compute conserved quantities
    system.conservation_entropy();
    system.conservation_BSQ();

    //===================================
    // print updated system and status
    io.print_conservation_status();
    io.print_system_state();
  }

  sw.Stop();
  formatted_output::summarize("All timesteps finished in "
                              + to_string(sw.printTime()) + " s");
}


// not yet defined
void BSQHydro::find_freeze_out_surface(){}


// not yet defined
void BSQHydro::print_results(){}
