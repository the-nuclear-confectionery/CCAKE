#include "BSQHydro.h"

using namespace ccake;

//Template instantiations
template class BSQHydro<1,EoM_default>;
template class BSQHydro<2,EoM_default>;
template class BSQHydro<3,EoM_default>;
template class BSQHydro<1,EoM_cartesian>;
template class BSQHydro<2,EoM_cartesian>;
template class BSQHydro<3,EoM_cartesian>;

/// @brief Constructor for the BSQHydro class.
/// @details This constructor initializes a BSQHydro object with the given
/// settings pointer. It also creates a system state object and initializes the
/// SPH workstation with the given settings and system state.
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
  wsPtr = std::make_shared<SPHWorkstation<D,TEOM>>(settingsPtr,systemPtr); 

  outPtr = std::make_shared<Output<D>>(settingsPtr,systemPtr);

  return;
}


/// @brief Shell function to call appropriate initial conditions readers
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
        if(e > settingsPtr->e_cutoff) systemPtr->add_particle( p );
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
/// @tparam D The number of spatial dimensions.
/// @tparam TEOM The equation of motion used in the simulation.
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
/// eight space-separated values. The first and fifth values are ignored by
/// ccake. These can be used by the IC generator to store additional
/// information. The second, third, and fourth values are the step sizes in the
/// x, y, and eta dimensions, respectively. The sixth, seventh, and eighth
/// values are the minimum values in the x, y, and eta dimensions, respectively.
/// A valid header would look like this:
///
/// #0 0.1 0.1 0.1 0 -10.0 -10.0 -0.5
///
/// The remaining lines are assumed to be in a space-separated format where each
/// line represents the properties of a single particle. The order of the
/// properties in each line should be as follows:
///
/// \f$x\f$ \f$y\f$ \f$\eta_s\f$ \f$\varepsilon\f$ \f$\rho_B\f$ \f$\rho_S\f$
/// \f$\rho_Q\f$ \f$u^x\f$ \f$u^y\f$ \f$u^{\eta}\f$ \f$\Pi\f$ \f$\pi^{xx}\f$
/// \f$\pi^{xy}\f$ \f$\pi^{x\eta}\f$ \f$\pi^{yy}\f$ \f$\pi^{y\eta}\f$
/// \f$\pi^{\eta\eta}\f$
///
/// where:
///
/// - \f$x\f$, \f$y\f$ and \f$\eta_s\f$ are the spatial coordinates of the
/// particle in the \f$x\f$, \f$y\f$ and \f$\eta_s\f$ directions, respectively.
/// - \f$\varepsilon\f$ is the energy density of the particle in units of 
/// GeV/fm\f$^3\f$.
/// - \f$\rho_B\f$, \f$\rho_S\f$ and \f$\rho_Q\f$ are the baryon, strangeness,
/// and electric charge densities of the particle, respectively. All densities 
/// are in units of 1/fm\f$^3\f$.
/// - \f$u^x\f$ \f$u^y\f$ \f$u^{\eta}\f$ are the fluid velocity components in 
/// the \f$x\f$, \f$y\f$ and \f$\eta_s\f$ directions, respectively.
/// - \f$\Pi\f$ is the bulk pressure of the fluid in units of GeV/fm\f$^3\f$.
/// - \f$\pi^{xx}\f$, \f$\pi^{xy}\f$, \f$\pi^{x\eta}\f$, \f$\pi^{yy}\f$,
/// \f$\pi^{y\eta}\f$, \f$\pi^{\eta\eta}\f$ are the components of the shear
/// tensor of the fluid.
///
/// @todo A suggestion is to add to the header additional info like
/// colliding system (AuAu, PbPb etc), colliding energy, impact 
/// parameter, etc. These would then be used as header of the outputs.
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
  if(settingsPtr->input_as_entropy==true){
    if(settingsPtr->e_cutoff > 0.0){
      std::cout<< "Using entropy as input, using ideal pion gas for e_cutoff" << std::endl;
    }
  }
  if (infile.is_open())
  {
    string line;
    int count_header_lines = 0;
    double x, y, eta, e, rhoB, rhoS, rhoQ, ux, uy, ueta, bulk, pixx, pixy, pixeta, piyy, piyeta, pietaeta;
    double ignore, stepX, stepY, stepEta, xmin, ymin, etamin;
    double max_e = 0.0;
    double min_e = 1e6;
    double s=0.0; //for input as entropy
    while (getline (infile, line))
    {
      istringstream iss(line);
      if(count_header_lines < total_header_lines)
      {
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
        if(settingsPtr->input_as_entropy==true){
          //convert the input to entropy density
          iss >> x >> y >> eta >> s >> rhoB >> rhoS >> rhoQ >> ux >> uy >> ueta >> bulk >> pixx >> pixy >> pixeta >> piyy >> piyeta >> pietaeta;
        }
        else{
          iss >> x >> y >> eta >> e >> rhoB >> rhoS >> rhoQ >> ux >> uy >> ueta >> bulk >> pixx >> pixy >> pixeta >> piyy >> piyeta >> pietaeta;
          e /= hbarc_GeVfm; // 1/fm^4
        }

        pixx /= hbarc_GeVfm; // 1/fm^4
        pixy /= hbarc_GeVfm; // 1/fm^4
        pixeta /= hbarc_GeVfm; // 1/fm^5
        piyy /= hbarc_GeVfm; // 1/fm^4
        piyeta /= hbarc_GeVfm; // 1/fm^5
        pietaeta /= hbarc_GeVfm; // 1/fm^6
        bulk /= hbarc_GeVfm; // 1/fm^4
        Particle<D> p;
        switch (D)
        {
          case 1:
            p.r(0) = eta;
            p.hydro.u(0) = ueta;
            p.hydro.shv(1,1) = pixx;
            p.hydro.shv(1,2) = pixy;
            p.hydro.shv(1,3) = pixeta;
            p.hydro.shv(2,1) = pixy;
            p.hydro.shv(3,1) = pixeta;
            p.hydro.shv(2,2) = piyy;
            p.hydro.shv(2,3) = piyeta;
            p.hydro.shv(3,2) = piyeta;
            p.hydro.shv(3,3) = pietaeta;
            break;
          case 2:
            p.r(0) = x;
            p.r(1) = y;
            p.hydro.u(0) = ux;
            p.hydro.u(1) = uy;
            p.hydro.shv(1,1) = pixx;
            p.hydro.shv(1,2) = pixy;
            p.hydro.shv(1,3) = pixeta;
            p.hydro.shv(2,1) = pixy;
            p.hydro.shv(3,1) = pixeta;
            p.hydro.shv(2,2) = piyy;
            p.hydro.shv(2,3) = piyeta;
            p.hydro.shv(3,2) = piyeta;
            p.hydro.shv(3,3) = pietaeta;
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
        if(settingsPtr->input_as_entropy==true){
          p.input.s = s;
        }
        else{
          p.input.e = e;
        }
        p.input.rhoB = rhoB;
        p.input.rhoS = rhoS;
        p.input.rhoQ = rhoQ;
        p.hydro.bulk = bulk;
        if(settingsPtr->input_as_entropy==true){
          
          //C = g_pi * pi^2 / (90) 
          double C = 3.0 * pow(M_PI,2.) / (90.0);
          double e_check = (3.*C)*pow(1./(4.*C),4./3.) * pow(s,4./3.);
          if(e_check > settingsPtr->e_cutoff) systemPtr->add_particle( p );
        }
        else{
          if(e > settingsPtr->e_cutoff) systemPtr->add_particle( p );
        }
        
        #ifdef DEBUG
        outfile << x << " " << y << " " << eta << " " << e*hbarc_GeVfm << " " << rhoB
                << " " << rhoS << " " << rhoQ << " " << ux << " " << uy << " " << ueta << endl;
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
/// @brief Shell function to initialize the hydrodynamics.
/// @details This function initializes the hydrodynamics by performing the 
/// following steps
/// - Initializing the workstation.
/// - Processing the initial conditions.
/// - Allocating memory for particles in the device.
/// - Finding the location of each particle in the phase diagram.
/// - Creating the linked list data structure.
/// - Setting up freeze-out procedure.
/// - Implementing initial smoothing required by the SPH formalism.
/// - Absorbing non-equilibrium pressure correction into bulk viscous pressure.
/// @see SPHWorkstation::initialize
/// @see SPHWorkstation::process_initial_conditions
/// @see SPHWorkstation::allocate_cabana_particles
/// @see SPHWorkstation::initialize_entropy_and_charge_densities
/// @see SPHWorkstation::initialize_linklist
/// @see SPHWorkstation::setup_freeze_out
/// @see SPHWorkstation::initial_smoothing
/// @see SPHWorkstation::set_bulk_Pi
/// @tparam D The number of spatial dimensions.
/// @tparam TEOM The equation of motion used in the simulation.
template<unsigned int D, template<unsigned int> typename TEOM>
void BSQHydro<D,TEOM>::initialize_hydrodynamics()
{
  formatted_output::announce("Initializing hydrodynamics");
  Stopwatch sw;
  sw.Start();

  // initialize workstation
  wsPtr->initialize();

  // trim initial conditions with low-energy density cut-off,
  // filling out initial conditions, and imposing initial freeze-out
  wsPtr->process_initial_conditions();


  // allocate memory for particles in device
  systemPtr->allocate_cabana_particles();
  

  // Create linked list data structure
  systemPtr->initialize_linklist();
  // for each particle, find location in phase diagram
  wsPtr->initialize_entropy_and_charge_densities();

  //Setup freeze-out
  wsPtr->setup_freeze_out();

  wsPtr->calculate_gamma_and_velocities();
  // implement initial smoothing required by SPH formalism
  wsPtr->initial_smoothing();

  // if initializing from full Tmunu, absorb non-equilibrium
  // pressure correction into bulk viscous pressure Pi
  wsPtr->set_bulk_Pi();
  //calculate extensive shear tensor
  wsPtr->calculate_extensive_shv();
  // calculate initial diffusion 
  wsPtr->set_diffusion();

  sw.Stop();
  formatted_output::report("hydrodynamics initialization finished in "
                              + to_string(sw.printTime()) + " s");
  

  #ifdef DEBUG_SLOW
  systemPtr->copy_device_to_host();
  std::ofstream thermo_file;
  thermo_file.open("initial_thermo.dat");
  thermo_file << "eta thermo.e(GeV/fm^3) input.e(GeV/fm^3) thermo.s(fm^-3) input.s(fm^-3)" << std::endl;
  for (auto & p : systemPtr->particles){
    //Print initial conditions
    for (int i = 0; i < D; i++) thermo_file << p.r(i) << " ";
    thermo_file << p.thermo.e*hbarc_GeVfm << " " << p.input.e*hbarc_GeVfm << " " << p.thermo.s << " " << p.input.s << " " << std::endl;
    //thermo_file << p.thermo.e*hbarc_GeVfm << " " << p.input.e*hbarc_GeVfm << " " << p.thermo.s << " " << p.thermo.rhoB << " " << p.thermo.rhoS 
    //            << " " << p.thermo.rhoQ << " ";
    //for (int i = 0; i < D; i++) thermo_file << p.hydro.u(i) << " ";
    //thermo_file << p.hydro.Bulk << " ";
    //for(int i = 0; i < D; i++){
    //  for(int j = i; j < D; j++){
    //    thermo_file << p.hydro.shv(i,j) << " ";
    //  }
    //}
    //thermo_file << p.hydro.shv33 << " " << std::endl;
  }
  thermo_file.close();
  exit(1);
  #endif

  return;
}

/// @brief Shell function to run the hydrodynamic evolution.
/// @details This function advances the hydrodynamic evolution by calling the
/// advance_timestep function of the workstation. It also computes the
/// conserved quantities, the eccentricities, and controls printing
/// of the system state, conservation status, and freeze-out surface.
/// @see SPHWorkstation::advance_timestep
/// @see output::print_system_state
/// @see output::print_freeze_out
/// @tparam D The number of spatial dimensions.
/// @tparam TEOM The equation of motion used in the simulation.
template<unsigned int D, template<unsigned int> typename TEOM>
void BSQHydro<D,TEOM>::run()
{
  formatted_output::announce("Beginning hydrodynamic evolution");
  Stopwatch sw;
  sw.Start();

  #ifdef DEBUG
  std::ofstream outfile;
  outfile.open("conservation.dat");
  #endif
  //===================================
  // initialize conserved quantities, etc.
  if (settingsPtr->print_conservation_status) {
    systemPtr->conservation_entropy(true);
    systemPtr->conservation_BSQ(true);
  }
  if (settingsPtr->calculate_observables) systemPtr->compute_eccentricities();


  //===================================
  // print initialized system and status
  if (settingsPtr->print_conservation_status) outPtr->print_conservation_status();
  #ifdef DEBUG
  outfile << systemPtr->t << " " << systemPtr->Eloss << " " << systemPtr->S << endl;
  #endif
  outPtr->print_system_state();

  //===================================
  // evolve until simulation terminates
  #ifdef DEBUG
  std::ofstream file;
  file.open("probe.dbg", std::ios::out | std::ios::trunc);
  if (file.is_open()) {
    file << "" ;
    file.close();
  }
  file.open("probe2.dbg", std::ios::out | std::ios::trunc);
  if (file.is_open()) {
    file << "" ;
    file.close();
  }
  #endif
  // Reset freeze-out file if performing freeze-out
  if (settingsPtr->particlization_enabled) {
    std::ofstream fo_file;
    fo_file.open(outPtr->get_freeze_out_filename(), std::ios::out | std::ios::trunc);
    if (fo_file.is_open()) {
      fo_file << "" ;
      fo_file.close();
    }
  }
  while ( wsPtr->continue_evolution() )
  {
    //===================================
    // workstation advances by given
    // timestep at given RK order
    wsPtr->advance_timestep( settingsPtr->dt, settingsPtr->rk_order );
    //wsPtr->regulator();

    //===================================
    // re-compute conserved quantities, etc.
    if (settingsPtr->print_conservation_status){
      systemPtr->conservation_entropy();
      systemPtr->conservation_BSQ();
    }
    if (settingsPtr->calculate_observables) systemPtr->compute_eccentricities();

    //===================================
    // print updated system and status
    if (settingsPtr->print_conservation_status) outPtr->print_conservation_status();
    #ifdef DEBUG
    outfile << systemPtr->t << " " << systemPtr->Eloss << " " << systemPtr->S << endl;
    #endif
    //outPtr->print_system_state();
    if (settingsPtr->hdf_evolution || settingsPtr->txt_evolution) 
    {
      if (systemPtr->number_of_elapsed_timesteps%100 == 0) outPtr->print_system_state();
    }
    if (settingsPtr->particlization_enabled) outPtr->print_freeze_out(wsPtr->freezePtr);

  }
  #ifdef DEBUG
  outfile.close();
  #endif

  sw.Stop();
  formatted_output::summarize("All timesteps finished in "
                              + to_string(sw.printTime()) + " s");
}
