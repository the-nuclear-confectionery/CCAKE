#include "BSQHydro.h"

using namespace ccake;

//Template instantiations
template class BSQHydro<1,EoM_default>;
template class BSQHydro<2,EoM_default>;
template class BSQHydro<3,EoM_default>;

/// @brief Constructor for the BSQHydro class.
/// @details This constructor initializes a BSQHydro object
/// with the given settings pointer. It also creates a
/// system state object and initializes the SPH workstation with the
/// given settings and system state.
/// @tparam D Dimensionality of the system.
/// @tparam TEOM Template for the Equation of Motion.
/// @param settingsPtr_in Shared pointer to the Settings object.
template<unsigned int D, template<unsigned int> typename TEOM>
BSQHydro<D,TEOM>::BSQHydro(std::shared_ptr<Settings> settingsPtr_in)
{
  //Initialize the settings pointer
  settingsPtr = settingsPtr_in;

  //Create the system state object - Necessary for reading in initial conditions
  systemPtr = std::make_shared<SystemState<D>>(settingsPtr);

  //Initialize the workstation
  wsPtr = std::make_shared<SPHWorkstation<D,TEOM>>(settingsPtr,systemPtr); //\TODO: If ever new EoM are implemented,
                                        // a switch case should be added here.
  
  outPtr = std::make_shared<Output<D>>(settingsPtr,systemPtr);

  return;
}



/// @brief Reads in initial conditions from a file
/// @details Reads in initial conditions from a file. The file format
/// is determined by the initial condition type specified in the
/// settings input. Currently supported types are:
///         - ICCING
///         - ccake
/// @tparam D TDimensionality of the system.
/// @tparam TEOM Template for the Equation of Motion.
template <unsigned int D, template <unsigned int> class TEOM>
void BSQHydro<D,TEOM>::read_in_initial_conditions(){

  formatted_output::report("Reading in initial conditions for hydrodynamics");

  string initial_condition_type = settingsPtr->IC_type;
  formatted_output::update("Initial conditions type: " + settingsPtr->IC_type);

  string IC_file = settingsPtr->IC_file;

  //----------------------------------------------------------------------------
  if (initial_condition_type == "ICCING")
  {
    read_ICCING();
  }
  else if (initial_condition_type == "ccake")
  {
    read_ccake();
  }

  return;
}

/// @brief Reads in the initial conditions for hydrodynamics from an ICCING file.
/// @details This function reads in the initial conditions for hydrodynamics
/// from an ICCING file. It extracts the particle position, energy density
/// and adds them to the system state object. The file path is specified in
/// the settings input file.
///
/// This is an specialization of the template function for 2D systems.
/// Currently, ICCING initial conditions are only available for 2D systems and
/// default EoM.
template<>
void BSQHydro<2,EoM_default>::read_ICCING()
{
  int total_header_lines = 1;
  string IC_file = settingsPtr->IC_file;

  ifstream infile(IC_file.c_str());
  formatted_output::update("Initial conditions file: " + IC_file);
  if (infile.is_open())
  {
    string line;
    int count_header_lines = 0;
    double x, y, e, rhoB, rhoS, rhoQ;
    double ignore, stepX, stepY, xmin, ymin;

    while (getline (infile, line))
    {
      istringstream iss(line);
      if(count_header_lines < total_header_lines)
      {
        settingsPtr->headers.push_back(line);
        iss >> ignore >> stepX >> stepY >> ignore >> xmin >> ymin;
        settingsPtr->stepx = stepX;
        settingsPtr->stepy = stepY;
        settingsPtr->xmin  = xmin;
        settingsPtr->ymin  = ymin;
        count_header_lines++;
      }
      else
      {
        iss >> x >> y >> e >> rhoB >> rhoS >> rhoQ;
        e /= hbarc_GeVfm;
        double ux = 0.0, uy = 0.0;

        Particle<2> p;
        p.r(0)       = x;
        p.r(1)       = y;
        p.input.e    = e;
        p.input.rhoB = rhoB;
        p.input.rhoS = rhoS;
        p.input.rhoQ = rhoQ;
        p.hydro.u(0) = ux;
        p.hydro.u(1) = uy;

        systemPtr->add_particle( p );
      }
    }
    infile.close();
  }
  else
  {
    std::cerr << "File " << IC_file << " could not be opened!\n";
    abort();
  }
}

/// @brief Reads in the initial conditions for hydrodynamics from an ICCING file.
/// @details This function should read in the initial conditions for hydrodynamics
/// from an ICCING file. It extracts the particle position, energy density
/// and adds them to the system state object. The file path is specified in
/// the settings input file.
///
/// This is the general template for reading ICCING files. Since ICCING
/// currently is 2D only, ww abort execution if the non-specialized function is,
/// called.
template<unsigned int D, template<unsigned int> typename TEOM>
void BSQHydro<D,TEOM>::read_ICCING()
{
    std::cerr << "ICCING initial conditions are only available for 2D simulations!\n";
    abort();
}


/// @brief Read the initial conditions from a file in the format of the CCAKE.
///
/// @details This function reads the initial conditions of a hydrodynamic
/// simulation from a file in the CCAKE format. The file format is as follows:
///
/// The first line is a header that contains information about the grid used to
/// generate the IC. It must start with a '#' character. It is expected to have
/// eight space-separated values. The first and fifth values are ignored, while
/// the second, third, and fourth values are the step sizes in the x, y, and eta
/// dimensions, respectively. The sixth, seventh, and eighth values are the
/// minimum values in the x, y, and eta dimensions, respectively. A valid header
/// would look like this:
///
/// #0 0.1 0.1 0.1 0 -10.0 -10.0 -0.5
///
/// The remaining lines are assumed to be in a space-separated format where each
/// line represents the properties of a single particle. The order of the
/// properties in each line should be as follows:
///
/// $x$ $y$ $\eta_s$ $\vareepsilon$ $\rho_B$ $\rho_S$ $\rho_Q$ $u^x$ $u^y$ $u^{\eta}$ $\Pi$ $\pi^{xx}$ $\pi^{xy}$ $\pi^{x\eta}$ $\pi^{yy}$ $\pi^{y\eta}$ $\pi^{\eta\eta}$
///
/// where:
///
/// - $x$, $y$ and $\eta_s$ are the spatial coordinates of the particle in the
///   $x$, $y$ and $\eta_s$ directions, respectively.
/// - $\vareepsilon$ is the energy density of the particle in units of GeV/fm$^3$.
/// - $\rho_B$, $\rho_S$ and $\rho_Q$ are the baryon, strangeness, and electric
///   charge densities of the particle, respectively. All densities are in units
///   of 1/fm$^3$.
/// - $u^x$ $u^y$ $u^{\eta}$ are the fluid velocity components in the $x$, $y$
///   and $\eta_s$ directions, respectively.
/// - $\Pi$ is the bulk pressure of the fluid in units of GeV/fm$^3$.
/// - $\pi^{xx}$, $\pi^{xy}$, $\pi^{x\eta}$, $\pi^{yy}$, $\pi^{y\eta}$,
///   $\pi^{\eta\eta}$ are the components of the shear tensor of the fluid.
///
/// @tparam D The number of spatial dimensions.
/// @tparam TEOM The equation of motion used in the simulation.
template<unsigned int D, template<unsigned int> typename TEOM>
void BSQHydro<D,TEOM>::read_ccake()
{
  int total_header_lines = 1;
  string IC_file = settingsPtr->IC_file;

  ifstream infile(IC_file.c_str());
  #ifdef DEBUG
  ofstream outfile;
  outfile.open("initial_conditions.dat");
  #endif
  formatted_output::update("Initial conditions file: " + IC_file);
  if (infile.is_open())
  {
    string line;
    int count_header_lines = 0;
    double x, y, eta, e, rhoB, rhoS, rhoQ, ux, uy, ueta, Bulk, pixx, pixy, pixeta, piyy, piyeta, pietaeta;
    double ignore, stepX, stepY, stepEta, xmin, ymin, etamin;

    while (getline (infile, line))
    {
      istringstream iss(line);
      if(count_header_lines < total_header_lines)
      {
        ///\todo Maybe we could store on the header additional info like
        ///      colliding system (AuAu, PbPb etc), colliding energy, impact parameter, etc.
        ///      These would then be used as header of the outputs.
        settingsPtr->headers.push_back(line);
        iss.ignore(256,'#');
        iss >> ignore >> stepX >> stepY >> stepEta >> ignore >> xmin >> ymin >> etamin;
        settingsPtr->stepx = stepX;
        settingsPtr->stepy = stepY;
        settingsPtr->stepEta = stepEta;
        settingsPtr->xmin  = xmin;
        settingsPtr->ymin  = ymin;
        settingsPtr->etamin = etamin;
        count_header_lines++;
      }
      else
      {
        iss >> x >> y >> eta >> e >> rhoB >> rhoS >> rhoQ >> ux >> uy >> ueta >> Bulk >> pixx >> pixy >> pixeta >> piyy >> piyeta >> pietaeta;
        e /= hbarc_GeVfm; // 1/fm^4
        pixx /= hbarc_GeVfm; // 1/fm^4
        pixy /= hbarc_GeVfm; // 1/fm^4
        pixeta /= hbarc_GeVfm; // 1/fm^5
        piyy /= hbarc_GeVfm; // 1/fm^4
        piyeta /= hbarc_GeVfm; // 1/fm^5
        pietaeta /= hbarc_GeVfm; // 1/fm^6
        Bulk /= hbarc_GeVfm; // 1/fm^4
        Particle<D> p;
        switch (D)
        {
          case 1:
            p.r(0) = eta;
            p.hydro.u(0) = ueta;
            p.hydro.shv(1,1) = pietaeta;
            break;
          case 2:
            p.r(0) = x;
            p.r(1) = y;
            p.hydro.u(0) = ux;
            p.hydro.u(1) = uy;
            p.hydro.shv(1,1) = pixx;
            p.hydro.shv(1,2) = pixy;
            p.hydro.shv(2,1) = pixy;
            p.hydro.shv(2,2) = piyy;
          case 3:
            p.r(0) = x;
            p.r(1) = y;
            p.r(2) = eta;
            p.hydro.u(0) = ux;
            p.hydro.u(1) = uy;
            p.hydro.u(2) = ueta;
            p.hydro.shv(1,1) = pixx;
            p.hydro.shv(1,2) = pixy;
            p.hydro.shv(1,3) = pixeta;
            p.hydro.shv(2,1) = pixy;
            p.hydro.shv(3,1) = pixeta;
            p.hydro.shv(2,2) = piyy;
            p.hydro.shv(2,3) = piyeta;
            p.hydro.shv(3,2) = piyeta;
            p.hydro.shv(3,3) = pietaeta;
        }
        p.hydro.shv33 = pietaeta;
        p.input.e    = e;
        p.input.rhoB = rhoB;
        p.input.rhoS = rhoS;
        p.input.rhoQ = rhoQ;
        p.hydro.Bulk = Bulk;
        systemPtr->add_particle( p );
        #ifdef DEBUG
        outfile << x << " " << y << " " << e*hbarc_GeVfm << " " << rhoB
                << " " << rhoS << " " << rhoQ << " " << ux << " " << uy << endl;
        #endif
      }
    }
    #ifdef DEBUG
    outfile.close();
    #endif
    infile.close();
  }
  else
  {
    std::cerr << "File " << IC_file << " could not be opened!" << std::endl;
    Kokkos::finalize();
    abort();
  }
}

////////////////////////////////////////////////////////////////////////////////
///
template<unsigned int D, template<unsigned int> typename TEOM>
void BSQHydro<D,TEOM>::initialize_hydrodynamics()
{
  formatted_output::announce("Initializing hydrodynamics");
  Stopwatch sw;
  sw.Start();

  // initialize workstation
  wsPtr->initialize();

  // initialize system state
  systemPtr->initialize();

  // trim initial conditions with low-energy density cut-off,
  // filling out initial conditions, and imposing initial freeze-out
  wsPtr->process_initial_conditions();

  // allocate memory for particles in device
  systemPtr->allocate_cabana_particles();

  // for each particle, find location in phase diagram
  wsPtr->initialize_entropy_and_charge_densities();

  // Create linked list data structure
  systemPtr->initialize_linklist();

  // implement initial smoothing required by SPH formalism
  wsPtr->initial_smoothing();

  // if initializing from full Tmunu, absorb non-equilibrium
  // pressure correction into bulk viscous pressure Pi
  wsPtr->set_bulk_Pi();

  sw.Stop();
  formatted_output::report("hydrodynamics initialization finished in "
                              + to_string(sw.printTime()) + " s");
  return;
}



////////////////////////////////////////////////////////////////////////////////
template<unsigned int D, template<unsigned int> typename TEOM>
void BSQHydro<D,TEOM>::run()
{
  formatted_output::announce("Beginning hydrodynamic evolution");
  Stopwatch sw;
  sw.Start();

  //===================================
  // initialize conserved quantities, etc.
  systemPtr->conservation_entropy(true);
  systemPtr->conservation_BSQ(true);
  systemPtr->compute_eccentricities();

  //===================================
  // print initialized system and status
  outPtr->print_conservation_status();
  outPtr->print_system_state();
  
  //===================================
  // evolve until simulation terminates
  int istep=0;
  while ( wsPtr->continue_evolution() )
  {
    //===================================
    // workstation advances by given
    // timestep at given RK order
    wsPtr->advance_timestep( settingsPtr->dt, rk_order );

    //===================================
    // re-compute conserved quantities, etc.
    //systemPtr->conservation_entropy();
    //systemPtr->conservation_BSQ();
    //systemPtr->compute_eccentricities();

    istep++;
    //===================================
    // print updated system and status
    //outPtr->print_conservation_status();
    if (istep%100 == 0) outPtr->print_system_state();
    if (istep==1000) break;
  }

  sw.Stop();
  formatted_output::summarize("All timesteps finished in "
                              + to_string(sw.printTime()) + " s");
  
}

/*
// not yet defined
template <unsigned int D>
void BSQHydro<D>::find_freeze_out_surface(){}


// not yet defined
template <unsigned int D>
void BSQHydro<D>::print_results(){}*/
