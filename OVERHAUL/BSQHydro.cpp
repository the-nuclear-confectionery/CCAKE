#include "BSQHydro.h"
#include "settings.h"
#include "input_output.h"
//#include "runge_kutta.h"
//#include "equations_of_motion.h"
#include "Stopwatch.h"

// Constructors and destructors.
BSQHydro::BSQHydro()
{

  // initialize I/O pointers
  io.set_EquationOfStatePtr( &eos );
  io.set_SettingsPtr( &settings );
  io.set_SystemStatePtr( &system );

  // initialize SPH workstation
  ws.set_EquationOfStatePtr( &eos );
  ws.set_SystemStatePtr( &system );
  ws.set_SettingsPtr( &settings );
  
  // initialize system state
  system.set_EquationOfStatePtr( &eos );
  system.set_SettingsPtr( &settings );

  // initialize EoS pointers
  eos.set_SettingsPtr( &settings );
  
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

  // initialize equation of state
  eos.init();

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

  return;
}


////////////////////////////////////////////////////////////////////////////////
void BSQHydro::run()
{
  cout << "Ready to start hydrodynamics\n";

  settings.t = settings.t0;

  // print the initial timestep
  io.print_system_state();
  cout << "Printed first timestep" << endl;

  system.conservation_entropy();
  system.conservation_BSQ();

  cout << setw(12) << setprecision(10)
       << "t=" << system.t << " S=" << system.S 
       << " " << system.Btotal << " " << system.Stotal
       << " " << system.Qtotal << endl;


  cout << "Now let's do the main evolution!" << endl;
  system.Ez = 0.0;

  while ( (system.t<settings.tend) && (system.number_part<system.n()) )
  {
    system.cfon = 1;

    // workstation advances by given timestep at given RK order
    ws.advance_timestep( settings.dt, rk_order );

    system.conservation_entropy();
    system.conservation_BSQ();

    // print energy/entropy and conserved charge totals
    cout << setw(12) << setprecision(10)
         << "t=" << system.t << " " << system.Eloss << " " << system.E0
         << " " << system.Etot << " " << system.S
         << " " << system.Btotal << " " << system.Stotal
         << " " << system.Qtotal << endl;

    // print system state, once per timestep
    io.print_system_state();

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // NOW DOING THIS AFTER ADDING FREEZE == 5 OPTION!!!!!!!!!!!!!!!!!
    system.number_part = system.get_frozen_out_count();
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    cout << "Check termination conditions: "
          << system.t << "   "
          << settings.tend << "   "
          << system.get_frozen_out_count() << "   "
          << system.number_part << "   "
          << system.n() << "   "
          << ( system.t < settings.tend ) << "   "
          << ( system.number_part < system.n() ) << endl;

  }
}


// not yet defined
void BSQHydro::find_freeze_out_surface(){}


// not yet defined
void BSQHydro::print_results(){}
