#include "BSQHydro.h"

namespace cc = ccake;

//Template instantiations
template class cc::BSQHydro<1>;
template class cc::BSQHydro<2>;
template class cc::BSQHydro<3>;

/// \brief Constructor for BSQHydro class
///
template<unsigned int D>
cc::BSQHydro<D>::BSQHydro(std::shared_ptr<Settings> settingsPtr_in)
{
  //Initialize the settings pointer
  settingsPtr = settingsPtr_in;


  // initialize I/O pointers
  io.set_EoS_type();
  io.set_EquationOfStatePtr( &ws.eos );
  io.set_SettingsPtr( settingsPtr.get() );
  io.set_SPHWorkstationPtr( &ws );
  io.set_SystemStatePtr( &system );

  // initialize SPH workstation
  ws.set_SystemStatePtr( &system );
  ws.set_SettingsPtr( settingsPtr.get() );

  // initialize system state
  system.set_SettingsPtr( settingsPtr.get() );

  return;
}


template <unsigned int D>
void cc::BSQHydro<D>::read_in_initial_conditions(){
  // tells Output to talk to system state and set initial system state
  io.read_in_initial_conditions();
  return;
}

////////////////////////////////////////////////////////////////////////////////
template <unsigned int D>
void cc::BSQHydro<D>::initialize_hydrodynamics()
{
  formatted_output::announce("Initializing hydrodynamics");
  Stopwatch sw;
  sw.Start();

  // initialize workstation
  ws.initialize();

  // initialize system state
  system.initialize();

  // trim initial conditions with low-energy density cut-off,
  // filling out initial conditions, and imposing initial freeze-out
  ws.process_initial_conditions();

  // sets nearest neighbors needed for smoothing
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
template <unsigned int D>
void cc::BSQHydro<D>::run()
{
  formatted_output::announce("Beginning hydrodynamic evolution");
  Stopwatch sw;
  sw.Start();

  //===================================
  // initialize conserved quantities, etc.
  system.conservation_entropy();
  system.conservation_BSQ();
  system.compute_eccentricities();

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
    ws.advance_timestep( settingsPtr->dt, rk_order );

    //===================================
    // re-compute conserved quantities, etc.
    system.conservation_entropy();
    system.conservation_BSQ();
    system.compute_eccentricities();

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
template <unsigned int D>
void cc::BSQHydro<D>::find_freeze_out_surface(){}


// not yet defined
template <unsigned int D>
void cc::BSQHydro<D>::print_results(){}
