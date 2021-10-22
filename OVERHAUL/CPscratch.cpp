



void BSQHydro::initialize_hydrodynamics()
{
  Imat.identity();  // need to set the identity matrix for use in simulation
  system.manualenter(ics,linklist);
  return;
}




void BSQHydro::run()
{
  // probably move this to initialization stuff
  system.initialize();
  /*for ( double t = system.tInit; t <= system.tInit; t += system.dt )
  {
    
  }*/
  system.BSQSimulation(); // dt and LinkList now members of system

  return;
}
