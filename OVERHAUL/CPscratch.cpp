



void BSQHydro::initialize_hydrodynamics()
{
  Imat.identity();  // need to set the identity matrix for use in simulation
  system.initialize();
  return;
}




void BSQHydro::run()
{
  system.BSQSimulation(); // dt and LinkList now members of system

  return;
}
