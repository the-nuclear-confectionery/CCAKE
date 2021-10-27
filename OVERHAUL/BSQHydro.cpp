#include "BSQHydro.h"
#include "settings.h"
#include "input_output.h"
#include "runge_kutta.h"
#include "equations_of_motion.h"
#include "Stopwatch.h"

// Constructors and destructors.
//BSQHydro::BSQHydro(){}
//BSQHydro::~BSQHydro(){}


void BSQHydro::load_settings_file( string path_to_settings_file )
{
  io.load_settings_file(path_to_settings_file); // sets the settings path in
  // InputOutput, then loads parameters into Input_parameters struct
  io.set_EoS_type(); // InputOutput talks to EoS and
  // tells it where to find its tables

  return;
}



void BSQHydro::set_results_directory( string path_to_results_directory )
{
  io.set_results_directory(path_to_results_directory); //set the results directory in InputOutput
  return;
}



void BSQHydro::read_in_initial_conditions()
{
  io.read_in_initial_conditions(); // tells InputOutput to talk to system state and set initial system state
  return;
}


void BSQHydro::trim_initial_conditions()
{
  vector<vector<double>> trimmed_grid;
  int cells_before_trim = initial_conditions.density_grid.size();
  for(int i=0; i<cells_before_trim; i++)
  {
    double e = initial_conditions.density_grid[i][2];
    double rhoB = initial_conditions.density_grid[i][3];
    double rhoS = initial_conditions.density_grid[i][4];
    double rhoQ = initial_conditions.density_grid[i][5];

    eos.s_out(e,rhoB,rhoS,rhoQ);
    if (eos.T() > input_parameters.Freeze_Out_Temperature)
    {
      trimmed_grid.push_back(initial_conditions.density_grid[i]);
    }
  }
  initial_conditions.density_grid = trimmed_grid;
}


////////////////////////////////////////////////////////////////////////////////
void BSQHydro::initialize_hydrodynamics()
{
  system.initialize();

  system.set_equation_of_state(eos);

  system.set_settings(settings);

  system.initialize_entropy_and_charge_densities(); // this should be in a switch/if

  system.initial_smoothing();

  return;
}


////////////////////////////////////////////////////////////////////////////////
void BSQHydro::run()
{
  cout << "Ready to start hydrodynamics\n";
  settings.frzc=0;
  settings.cf=0;

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

    cout << "t=" << system.t << " S=" << system.S 
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
    settings.cfon=1;


    cout << "Entering here:" << endl;

    RK::bsq_second_order( settings.dt, eom, system );
    system.conservation_entropy();
    system.conservation_BSQ();

    cout << "t=" << system.t << " " <<  system.Eloss << " " << system.S
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
        cout << "t=" << system.t << " S=" << system.S << endl;  // outputs time step
        //out.bsqsveprofile(system);   // energy density profile
        cout << "eloss= " << system.t << " " <<  system.Eloss << endl;
        //out.conservation(system); // conservation of energy
      }
    }

  }
}


void BSQHydro::find_freeze_out_surface(){}


void BSQHydro::print_results(){}
