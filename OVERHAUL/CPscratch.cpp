////////////////////////////////////////////////////////////////////////////////
void BSQHydro::initialize_equation_of_state()
{
  //////////////////////////////////////////////////////////////////////////////
  // SET EQUATION OF STATE
  // rewrite by C. Plumberg: allow for different EOS format if using BSQ
  double efcheck = 0.0, sfcheck = 0.0;
  if ( linklist.visc == 4 )	//if we're running BSQ (table is only option)
  {
    bool using_HDF = false;
    if (using_HDF)
    {
      string quantityFile   = ifolder + std::string("quantityFile.h5");
      string derivativeFile = ifolder + std::string("derivFile.h5");
      std::cout << "Using BSQ Equation of State table from: "
                << quantityFile << " and " << derivativeFile << "\n";

      eos.init( quantityFile, derivativeFile );
    }
    else
    {
      string quantityFilename   = "EoS_Taylor_AllMu_T0_1200.dat";
      string derivativeFilename = "EoS_Taylor_AllMu_Derivatives_T0_1200.dat";
      string quantityFile       = ifolder + quantityFilename;
      string derivativeFile     = ifolder + derivativeFilename;
      std::cout << "Using BSQ Equation of State table from: "
                << quantityFile << " and " << derivativeFile << "\n";

      eos.init( quantityFile, derivativeFile );
    }

    eos.eosin( eostype );			// does nothing!
    const double freeze_out_T_at_mu_eq_0
                  = 0.15/hbarc_GeVfm;	//1/fm
    efcheck       = EOS0.efreeze( freeze_out_T_at_mu_eq_0 );
    sfcheck       = EOS0.sfreeze( freeze_out_T_at_mu_eq_0 );
    //efcheck = 0.266112/0.1973;
    //sfcheck = 2.05743;

    std::cout << "efcheck = " << efcheck*hbarc_GeVfm << " GeV/fm^3\n";
    std::cout << "sfcheck = " << sfcheck << " 1/fm^3\n";
  }
  else
  {
    std::cerr << "This EoS model not currently supported!" << std::endl;
  }

  // Make sure Particle and SystemState classes can access equation of state
  Particle::set_equation_of_state(    &EOS0 );
  SystemState::set_equation_of_state( &EOS0 );


  return;
}