#include "BSQHydro.h"
#include "settings.h"
#include "input_output.h"
#include "runge_kutta.h"
#include "equations_of_motion.h"
#include "Stopwatch.h"

// Constructors and destructors.
BSQHydro::BSQHydro()
{

  // initialize I/O pointers
  io.set_EquationOfStatePtr( &eos );
  io.set_EquationsOfMotionPtr( &eom );
  io.set_SettingsPtr( &settings );
  io.set_SystemStatePtr( &system );

  // initialize SPH workstation
  ws.set_EquationOfStatePtr( &eos );
  ws.set_SystemStatePtr( &system );
  ws.set_SettingsPtr( &settings );
  
  // initialize system state
  system.set_EquationOfStatePtr( &eos );
  system.set_SettingsPtr( &settings );
  
  return;
}

//BSQHydro::~BSQHydro(){}


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
  //settings.frzc=0;
  //settings.cf=0;

  /*BBMG<2> bbmg(system);
  bbmg.initial(system);
  cout << "started bbmg" << endl;*/

  settings.t = settings.t0;

  if ( settings.qmf == 1 || settings.qmf == 3 )
  {
    //out.bsqsveprofile(system);
    cout << "printed first timestep" << endl;

    system.conservation_entropy();
    system.conservation_BSQ();

    cout << setw(12) << setprecision(10)
         << "t=" << system.t << " S=" << system.S 
         << " " << system.Btotal << " " << system.Stotal
         << " " << system.Qtotal << endl;

    if (settings.qmf==1) exit(0);
  }
  else if(settings.qmf==4)
  {
    //out.eccout(system);
    cout << "eccentricity printed" << endl;
    exit(0);
  }


  cout << "Now let's do the main evolution!" << endl;
  system.Ez=0;

  while ((system.t<settings.tend)&&(system.number_part<system.n()))
  {
    system.cfon = 1;

cout << "TEST LOOP: " << system.t << "   " << settings.tend
      << "   " << system.number_part << "   " << system.n() << endl;


    cout << "Entering here:" << endl;

    RK::bsq_second_order( settings.dt, eom, system, ws );
    system.conservation_entropy();
    system.conservation_BSQ();

    cout << setw(12) << setprecision(10)
         << "t=" << system.t << " " <<  system.Eloss << " " << system.S
         << " " << system.Btotal << " " << system.Stotal
         << " " << system.Qtotal <<  endl;

    //out.bsqsveprofile(system);


    //if (settings.cf>0) out.bsqsvFOprint(system);

    if (settings.qmf==3)
    {
      double tsub=system.t-floor(system.t);
      // if you add more points to print, must also change system<D>::setup and multiply steps=floor(tend-t0)+1; by the extra number of print offs / 1fm/c
      if (tsub<(0.0+system.dt*0.99)||(tsub>=1-+system.dt*0.99)) // uncomment if you want to observe energydensity profile, conservation of energy or do a Gubser check
      {
        system.conservation_entropy();
        cout << setw(12) << setprecision(10)
         << "t=" << system.t << " S=" << system.S << endl;  // outputs time step
        //out.bsqsveprofile(system);   // energy density profile
        cout << "eloss= " << system.t << " " <<  system.Eloss << endl;
        //out.conservation(system); // conservation of energy
      }
    }

    // print system state, once per timestep
    print_system_state();

  }
}


// not yet defined
void BSQHydro::find_freeze_out_surface(){}


// not yet defined
void BSQHydro::print_results(){}
