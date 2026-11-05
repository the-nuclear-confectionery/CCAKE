#include <iostream>

#include "BSQHydro.h"
#include "phase_profiler.h"
#include "wvt_glass.h"

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

  outPtr = std::make_shared<Output<D,TEOM>>(settingsPtr,systemPtr,wsPtr);

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
  if (settingsPtr->wvt_enabled)
  {
    read_wvt_glass();
  }
  else if (initial_condition_type == "ICCING")
  {
    read_ICCING();
  }
  else if (initial_condition_type == "ccake")
  {
    read_ccake();
  }
  else if (initial_condition_type == "freezein")
  {
    read_freezein();
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Builds an adaptive WVT glass from the gridded IC file and stages its rows
/// as particles through the normal "ccake" / "freezein" row conversion.
/// @details Replaces the former standalone wvt_ic_generator stage (see
/// include/wvt_glass.h). The glass is made before the cutoff, over the whole source;
/// energy_cutoff is applied afterwards by add_ccake_particle / add_freezein_particle,
/// as for a grid IC. The 3D stretch is eta_h / eta_h_eta with the factorized kernel
/// and 1 otherwise (unless wvt.eta_stretch is set), the anisotropy the kernel needs.
/// The grid header is set to the one the generator wrote for the glass (mean
/// spacing and glass minima); the neighbour grid bounds are built from it.
template <unsigned int D, template <unsigned int> class TEOM>
void BSQHydro<D,TEOM>::read_wvt_glass()
{
  const string& type = settingsPtr->IC_type;
  if (type != "ccake" && type != "freezein")
  {
    std::cerr << "ERROR: initial_conditions.wvt needs type ccake or freezein, not "
              << type << std::endl;
    Kokkos::finalize();
    abort();
  }

  wvt::GlassSettings glass_settings = settingsPtr->wvt_glass;
  if (glass_settings.eta_stretch <= 0.0)
    glass_settings.eta_stretch = (D == 3 && settingsPtr->factorized_kernel && settingsPtr->eta_h_eta > 0.0)
                               ? settingsPtr->eta_h / settingsPtr->eta_h_eta : 1.0;
  if (settingsPtr->smoothing_length_mode == "constant")
    formatted_output::report("WARNING: a WVT glass is density-adaptive, but "
                             "smoothing_length_mode is constant");

  const string source_file = (type == "freezein") ? settingsPtr->freezein_input_file.string()
                                                  : settingsPtr->IC_file.string();
  wvt::Glass glass = wvt::make_glass(source_file, glass_settings);
  if (glass.dimension != static_cast<int>(D))
  {
    std::cerr << "ERROR: the WVT source is " << glass.dimension << "D but the run is "
              << D << "D" << std::endl;
    Kokkos::finalize();
    abort();
  }
  const int needed_columns = (type == "freezein") ? 25
                           : (settingsPtr->input_initial_diffusion ? 29 : 17);
  if (glass.number_of_columns < needed_columns)
  {
    std::cerr << "ERROR: the WVT source has " << glass.number_of_columns
              << " columns, the " << type << " reader needs " << needed_columns << std::endl;
    Kokkos::finalize();
    abort();
  }
  if (!settingsPtr->wvt_output_file.empty())
    wvt::write_glass(settingsPtr->wvt_output_file, glass);

  settingsPtr->stepx   = glass.step_x;
  settingsPtr->stepy   = glass.step_y;
  settingsPtr->stepEta = glass.step_eta;
  settingsPtr->xmin    = glass.x_min;
  settingsPtr->ymin    = glass.y_min;
  settingsPtr->etamin  = glass.eta_min;

  std::size_t n_kept = 0;
  for (int i = 0; i < glass.number_of_rows(); ++i)
  {
    if (type == "freezein") { if (add_freezein_particle(glass.row(i))) ++n_kept; }
    else                    { add_ccake_particle(glass.row(i)); }
  }
  if (type == "freezein")
    formatted_output::update("Freeze-in: kept " + std::to_string(n_kept) + " of "
                             + std::to_string(glass.number_of_rows())
                             + " glass particles above e_cutoff = "
                             + std::to_string(settingsPtr->e_cutoff) + " GeV/fm^3");
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
        (settingsPtr->baryon_charge_enabled) ? p.input.rhoB = rhoB : p.input.rhoB = 0.0;
        (settingsPtr->strange_charge_enabled) ? p.input.rhoS = rhoS : p.input.rhoS = 0.0;
        (settingsPtr->electric_charge_enabled) ? p.input.rhoQ = rhoQ : p.input.rhoQ = 0.0;        
        if(e > settingsPtr->e_cutoff/hbarc_GeVfm) systemPtr->add_particle( p );
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
    // One data row (layout in add_ccake_particle). Declared outside the loop, like
    // the named locals it replaces, so a short row keeps the previous row's values.
    double row[29] = {0.0};
    const int number_of_columns = settingsPtr->input_initial_diffusion ? 29 : 17;
    double ignore, stepX, stepY, stepEta, xmin, ymin, etamin;
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
        // Skip blank lines and any extra comment lines in the data section.
        // Some ICs (e.g. SMASH-derived ic_hy files) carry several '#' metadata
        // lines after the single grid header; without this, those lines would be
        // parsed as data, leaving uninitialized x/y/eta/e -> NaN/off-grid particles.
        const size_t nb = line.find_first_not_of(" \t\r\n");
        if (nb == std::string::npos || line[nb] == '#') continue;

        for (int c = 0; c < number_of_columns; ++c) iss >> row[c];
        add_ccake_particle(row);

        #ifdef DEBUG
        for (int c = 0; c < 10; ++c) outfile << row[c] << (c < 9 ? " " : "\n");
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
/// @brief Converts one row of a "ccake" IC and adds it as a particle if it is
/// above the energy cutoff.
/// @details Row layout: x y eta e(or s) rhoB rhoS rhoQ ux uy ueta bulk pixx pixy
/// pixeta piyy piyeta pietaeta, then qB0..3 qS0..3 qQ0..3 when
/// input_initial_diffusion is set. Energies and stresses are in GeV/fm^3.
/// Used for rows read from the file and for rows of a WVT glass.
template<unsigned int D, template<unsigned int> typename TEOM>
void BSQHydro<D,TEOM>::add_ccake_particle(const double* row)
{
  const double x = row[0], y = row[1], eta = row[2];
  double e = 0.0, s = 0.0;
  if(settingsPtr->input_as_entropy==true){
    s = row[3];
  }
  else{
    e = row[3] / hbarc_GeVfm; // 1/fm^4
  }
  const double rhoB = row[4], rhoS = row[5], rhoQ = row[6];
  const double ux = row[7], uy = row[8], ueta = row[9];
  const double bulk     = row[10] / hbarc_GeVfm; // 1/fm^4
  const double pixx     = row[11] / hbarc_GeVfm; // 1/fm^4
  const double pixy     = row[12] / hbarc_GeVfm; // 1/fm^4
  const double pixeta   = row[13] / hbarc_GeVfm; // 1/fm^5
  const double piyy     = row[14] / hbarc_GeVfm; // 1/fm^4
  const double piyeta   = row[15] / hbarc_GeVfm; // 1/fm^5
  const double pietaeta = row[16] / hbarc_GeVfm; // 1/fm^6
  const double qB0 = row[17], qB1 = row[18], qB2 = row[19], qB3 = row[20];
  const double qS0 = row[21], qS1 = row[22], qS2 = row[23], qS3 = row[24];
  const double qQ0 = row[25], qQ1 = row[26], qQ2 = row[27], qQ3 = row[28];

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
      break;
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
  (settingsPtr->baryon_charge_enabled) ? p.input.rhoB = rhoB : p.input.rhoB = 0.0;
  (settingsPtr->strange_charge_enabled) ? p.input.rhoS = rhoS : p.input.rhoS = 0.0;
  (settingsPtr->electric_charge_enabled) ? p.input.rhoQ = rhoQ : p.input.rhoQ = 0.0;
  p.hydro.bulk = bulk;
  if(settingsPtr->input_initial_diffusion){
    if(settingsPtr->baryon_charge_enabled){
      p.hydro.diffusion(0,0) = qB0;
      p.hydro.diffusion(0,1) = qB1;
      p.hydro.diffusion(0,2) = qB2;
      p.hydro.diffusion(0,3) = qB3;
    } else {
      p.hydro.diffusion(0,0) = 0.0;
      p.hydro.diffusion(0,1) = 0.0;
      p.hydro.diffusion(0,2) = 0.0;
      p.hydro.diffusion(0,3) = 0.0;
    }
    if(settingsPtr->strange_charge_enabled){
      p.hydro.diffusion(1,0) = qS0;
      p.hydro.diffusion(1,1) = qS1;
      p.hydro.diffusion(1,2) = qS2;
      p.hydro.diffusion(1,3) = qS3;
    } else {
      p.hydro.diffusion(1,0) = 0.0;
      p.hydro.diffusion(1,1) = 0.0;
      p.hydro.diffusion(1,2) = 0.0;
      p.hydro.diffusion(1,3) = 0.0;
    }
    if(settingsPtr->electric_charge_enabled){
      p.hydro.diffusion(2,0) = qQ0;
      p.hydro.diffusion(2,1) = qQ1;
      p.hydro.diffusion(2,2) = qQ2;
      p.hydro.diffusion(2,3) = qQ3;
    } else {
      p.hydro.diffusion(2,0) = 0.0;
      p.hydro.diffusion(2,1) = 0.0;
      p.hydro.diffusion(2,2) = 0.0;
      p.hydro.diffusion(2,3) = 0.0;
    }
  }
  else{
    //set diffusion to zero
    p.hydro.diffusion(0,0) = 0.0;
    p.hydro.diffusion(0,1) = 0.0;
    p.hydro.diffusion(0,2) = 0.0;
    p.hydro.diffusion(0,3) = 0.0;
    p.hydro.diffusion(1,0) = 0.0;
    p.hydro.diffusion(1,1) = 0.0;
    p.hydro.diffusion(1,2) = 0.0;
    p.hydro.diffusion(1,3) = 0.0;
    p.hydro.diffusion(2,0) = 0.0;
    p.hydro.diffusion(2,1) = 0.0;
    p.hydro.diffusion(2,2) = 0.0;
    p.hydro.diffusion(2,3) = 0.0;
  }


  if(settingsPtr->input_as_entropy==true){
    
    //C = g_pi * pi^2 / (90) 
    double C = 3.0 * pow(M_PI,2.) / (90.0);
    double e_check = (3.*C)*pow(1./(4.*C),4./3.) * pow(s,4./3.);
    if(e_check > settingsPtr->e_cutoff/hbarc_GeVfm) systemPtr->add_particle( p );
  }
  else{
    if(e > settingsPtr->e_cutoff/hbarc_GeVfm) systemPtr->add_particle( p );
  }
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Reads a raw (un-diagonalized) contravariant T^{mu nu} grid produced by
/// AMPTGenesis (output.raw_tmunu) and stages each above-cutoff cell as an SPH
/// particle for the energy-momentum-conserving freeze-in matching.
/// @details The smearing of AMPT partons into T^{mu nu} (+ lab-frame B/S/Q
/// 4-currents) is done by AMPTGenesis; this reader only ingests the grid and
/// stores the conserved surface densities T^{tau mu} and j_a^mu on each
/// particle. The actual matching to (epsilon, u^mu, rho_{B,S,Q}) — using the
/// EoS to close P = P(epsilon, n) — is performed later in
/// SPHWorkstation::freeze_in_match(), once the EoS is initialized.
///
/// File format (see AMPTGenesis.cpp write_raw_tmunu):
///   header : #0 dx dy deta 0 xmin ymin etamin
///   comment: # raw (un-diagonalized) contravariant T^{mu nu} ...
///   comment: # columns: x y eta Ttt Ttx Tty Ttn Txx Txy Txn Tyy Tyn Tnn
///                       jB0..3 jS0..3 jQ0..3
///   rows   : x y eta <10 T components> <12 current components>
/// Indices are 0=tau, 1=x, 2=y, 3=eta; T^{mu nu} is in GeV/fm^3 and the
/// currents in 1/fm^3.
/// @tparam D Dimensionality of the system.
/// @tparam TEOM Template for the Equation of Motion.
template <unsigned int D, template <unsigned int> class TEOM>
void BSQHydro<D,TEOM>::read_freezein()
{
  const string IC_file = settingsPtr->freezein_input_file.string();
  ifstream infile(IC_file.c_str());
  formatted_output::update("Freeze-in (raw T^{mu nu}) file: " + IC_file);

  if (!infile.is_open())
  {
    std::cerr << "Freeze-in file " << IC_file << " could not be opened!" << std::endl;
    Kokkos::finalize();
    abort();
  }

  string line;
  bool header_read = false;
  std::size_t n_cells = 0, n_kept = 0;
  // e_cutoff is stored in GeV/fm^3; T^{tau tau} is also in GeV/fm^3 here.
  const double e_cutoff = settingsPtr->e_cutoff;

  while (getline(infile, line))
  {
    // Skip blank lines.
    if (line.find_first_not_of(" \t\r\n") == string::npos) continue;

    // First '#'-line is the grid header; later '#'-lines are comments/counts.
    if (line[line.find_first_not_of(" \t")] == '#')
    {
      if (!header_read)
      {
        settingsPtr->headers.push_back(line);
        istringstream iss(line);
        double ignore, stepX, stepY, stepEta, xmin, ymin, etamin;
        iss.ignore(256, '#');
        iss >> ignore >> stepX >> stepY >> stepEta >> ignore >> xmin >> ymin >> etamin;
        settingsPtr->stepx   = stepX;
        settingsPtr->stepy   = stepY;
        settingsPtr->stepEta = stepEta;
        settingsPtr->xmin    = xmin;
        settingsPtr->ymin    = ymin;
        settingsPtr->etamin  = etamin;
        header_read = true;
      }
      continue;
    }

    istringstream iss(line);
    double row[25];
    bool complete_row = true;
    for (int c = 0; c < 25 && complete_row; ++c) complete_row = static_cast<bool>(iss >> row[c]);
    if (!complete_row) continue; // not a data row (e.g. trailing comment we did not catch)

    ++n_cells;
    if (add_freezein_particle(row)) ++n_kept;
  }
  infile.close();

  formatted_output::update("Freeze-in: kept " + std::to_string(n_kept)
                           + " of " + std::to_string(n_cells)
                           + " cells above e_cutoff = "
                           + std::to_string(e_cutoff) + " GeV/fm^3");
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Stages one freeze-in row as a particle unless it is safely below the cutoff.
/// @details Row layout: x y eta Ttt Ttx Tty Ttn Txx Txy Txn Tyy Tyn Tnn jB0..3 jS0..3
/// jQ0..3 (T in GeV/fm^3, currents in 1/fm^3). Since the rest-frame energy density
/// satisfies epsilon <= T^{tau tau}, any row with T^{tau tau} <= e_cutoff is below the
/// cutoff; the decisive cut on epsilon follows in SPHWorkstation::freeze_in_match().
/// Used for rows read from the file and for rows of a WVT glass.
/// @return true if the particle was added.
template <unsigned int D, template <unsigned int> class TEOM>
bool BSQHydro<D,TEOM>::add_freezein_particle(const double* row)
{
  const double x = row[0], y = row[1], eta = row[2];
  const double Ttt = row[3], Ttx = row[4], Tty = row[5], Ttn = row[6];
  const double jB0 = row[13], jB1 = row[14], jB2 = row[15], jB3 = row[16];
  const double jS0 = row[17], jS1 = row[18], jS2 = row[19], jS3 = row[20];
  const double jQ0 = row[21], jQ1 = row[22], jQ2 = row[23], jQ3 = row[24];
  if (Ttt <= settingsPtr->e_cutoff) return false; // GeV/fm^3

  Particle<D> p;
  switch (D)
  {
    case 1:
      p.r(0) = eta;
      break;
    case 2:
      p.r(0) = x;
      p.r(1) = y;
      break;
    case 3:
      p.r(0) = x;
      p.r(1) = y;
      p.r(2) = eta;
  }

  // Stash the conserved surface densities (Milne, 0=tau,1=x,2=y,3=eta).
  p.freezein_Ttau[0] = Ttt;  p.freezein_Ttau[1] = Ttx;
  p.freezein_Ttau[2] = Tty;  p.freezein_Ttau[3] = Ttn;
  p.freezein_jB[0] = jB0; p.freezein_jB[1] = jB1; p.freezein_jB[2] = jB2; p.freezein_jB[3] = jB3;
  p.freezein_jS[0] = jS0; p.freezein_jS[1] = jS1; p.freezein_jS[2] = jS2; p.freezein_jS[3] = jS3;
  p.freezein_jQ[0] = jQ0; p.freezein_jQ[1] = jQ1; p.freezein_jQ[2] = jQ2; p.freezein_jQ[3] = jQ3;

  systemPtr->add_particle(p);
  return true;
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

  //if using source terms set them up
  if(settingsPtr->source_type != "disabled"){
    wsPtr->setup_source_terms();
  }
  
  wsPtr->calculate_gamma_and_velocities();
  // implement initial smoothing required by SPH formalism
  wsPtr->initial_smoothing();

  // initialize jets from particle temperatures and densities
  if(settingsPtr->jets_type != "disabled"){
    wsPtr->initialize_jets_bbmg();
  }
  // if initializing from full Tmunu, absorb non-equilibrium
  // pressure correction into bulk viscous pressure Pi
  //wsPtr->set_bulk_Pi();
  ////calculate extensive shear tensor
  ////wsPtr->calculate_extensive_shv();
  //// calculate initial diffusion 
  //wsPtr->set_diffusion();

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
  size_t n_steps_ran = 0;   // counts evolution timesteps (for steps/s under CCAKE_PROFILE)

  #ifdef DEBUG
  std::ofstream outfile;
  outfile.open("conservation.dat");
  #endif

  //====================================================================================
  // initialize conserved quantities, etc.
  // print conservation status to a file labelled "conservation.dat" in the output dir
  string out_dir = settingsPtr->results_directory;
  // std::ofstream outfile;
  if (settingsPtr->print_conservation_status) {
    systemPtr->conservation_entropy(true);
    systemPtr->conservation_BSQ(true);
    outPtr->print_conservation_status();
    outPtr->buffer_conservation_line();   // buffer initial state
  }

  //====================================================================================
  // printing eccentricities for given rapidity slices to respective files
  if (settingsPtr->calculate_observables) systemPtr->compute_eccentricities();

  // TODO:: Print to a file somehow at each time step!
  
  if (settingsPtr->hdf_evolution || settingsPtr->txt_evolution) {
    outPtr->print_system_state();
  }

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
  //====================================================================================
  // Printing neighbors for max and min particles based on entropy density
  // std::ofstream outfile_neighbors;
  std::vector<int> idx;
  // std::string n_path = out_dir + "/neighbors_" + std::to_string(idx) + ".dat";
  if(settingsPtr->get_neighbors){
    idx = systemPtr->entropy_density_based_idx_search();
    std::vector<std::string> labels = {"min", "max"};
    for (int i = 0; i < 2; ++i) {
        std::string n_path = out_dir + "/neighbors_" + labels[i] + ".dat";

        std::ofstream outfile_neighbors(n_path.c_str());
        outfile_neighbors << "t #neighbors x y eta" << std::endl;

        std::vector<std::array<double, 4>> result = systemPtr->get_particle_data(idx[i]);
        outfile_neighbors << systemPtr->t << " "
                          << result[0][0] << " "
                          << result[0][1] << " "
                          << result[0][2] << " "
                          << result[0][3] << std::endl;
    }
  }

  // optional benchmark cap: stop after CCAKE_MAX_STEPS evolution steps (0 = unlimited)
  const char* max_steps_env = std::getenv("CCAKE_MAX_STEPS");
  const size_t max_steps = (max_steps_env && std::atol(max_steps_env) > 0)
                           ? static_cast<size_t>(std::atol(max_steps_env)) : 0;
  const auto t_loop0 = std::chrono::high_resolution_clock::now();
  while ( wsPtr->continue_evolution() )
  {
    //===================================
    // workstation advances by given
    // timestep at given RK order
    wsPtr->advance_timestep( settingsPtr->dt, settingsPtr->rk_order );
    ++n_steps_ran;

    // interim per-phase report every CCAKE_PROFILE_EVERY steps (0 = off / end-only)
    {
      const int every = ccake::PhaseProfiler::instance().report_every();
      if ( every > 0 && (n_steps_ran % static_cast<size_t>(every)) == 0 )
      {
        const double el = std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - t_loop0).count();
        formatted_output::summarize("[profile] step " + to_string(n_steps_ran)
            + "   (" + to_string(el > 0.0 ? n_steps_ran/el : 0.0) + " steps/s running)");
        ccake::PhaseProfiler::instance().report(std::cout);
      }
    }

    if ( max_steps > 0 && n_steps_ran >= max_steps ) break;

    //==================================================================
    // re-compute conserved quantities, etc.
    // print updated system and status
    if (settingsPtr->print_conservation_status){
      systemPtr->conservation_entropy();
      systemPtr->conservation_BSQ();
      outPtr->print_conservation_status();
      outPtr->buffer_conservation_line();
    }
    if (settingsPtr->calculate_observables) systemPtr->compute_eccentricities();

    #ifdef DEBUG
    outfile << systemPtr->t << " " << systemPtr->Eloss << " " << systemPtr->S << endl;
    #endif

    // per-step causality summary (A, cheap) + strided minimal per-particle dump (B)
    if (settingsPtr->check_causality)
    {
      systemPtr->copy_device_to_host();
      outPtr->buffer_causality_line();
      if (settingsPtr->causality_minimal && (n_steps_ran % settingsPtr->causality_minimal_stride == 0))
        outPtr->print_causality_minimal_to_txt();
    }
    // full evolution output, optionally strided via evolution_stride
    if ( (settingsPtr->hdf_evolution || settingsPtr->txt_evolution || settingsPtr->jet_evolution)
         && (n_steps_ran % settingsPtr->evolution_stride == 0) )
    {
      outPtr->print_system_state();
    }
    if (settingsPtr->particlization_enabled) outPtr->print_freeze_out(wsPtr->freezePtr, systemPtr->Btotal0, systemPtr->Qtotal0, systemPtr->Stotal0);

    // =================================================================
    // Printing neighbors for each timestep
    if(settingsPtr->get_neighbors){
      std::vector<std::string> labels = {"min", "max"};
      for (int i = 0; i < 2; ++i) {
          std::string n_path = out_dir + "/neighbors_" + labels[i] + ".dat";
          std::ofstream outfile_neighbors(n_path.c_str(), std::ios::app);
          std::vector<std::array<double, 4>> result = systemPtr->get_particle_data(idx[i]);
          outfile_neighbors << systemPtr->t << " "
                            << result[0][0] << " "
                            << result[0][1] << " "
                            << result[0][2] << " "
                            << result[0][3] << std::endl;
      }
    }
  }
  // if(settingsPtr->get_neighbors){
  //   outfile_neighbors.close();
  // }
  // if (settingPtr->print_conservation_status) outfile.close();
  // #ifdef DEBUG
  // outfile.close();
  // #endif

  //=======================================================================
  // print observables for each timestep
  if (settingsPtr->calculate_observables) {
    std::ofstream outfile;
    for (int j = 0; j < systemPtr->eta_slices.size(); ++j){
      std::ostringstream eta_stream;
      eta_stream << std::fixed << std::setprecision(1) << systemPtr->eta_slices[j];
      string ecc_path = out_dir + "/eccentricities_" + eta_stream.str() + ".dat";
      outfile.open(ecc_path.c_str());
      outfile << "t " << "e_2_X " << "e_2_P " << "count_X " << "count_P " << endl;
      for (int i = 0; i < systemPtr->timesteps.size(); ++i)
      {
        outfile << systemPtr->timesteps[i] << " " \
        << systemPtr->e_2_X_history_by_slice[j][i] << " " \
        << systemPtr->e_2_P_history_by_slice[j][i] << " " \
        << systemPtr->count_X_history_by_slice[j][i] << " " \
        << systemPtr->count_P_history_by_slice[j][i] << endl;
      }
      outfile.close();
    }
  }

  sw.Stop();
  formatted_output::summarize("All timesteps finished in "
                              + to_string(sw.printTime()) + " s");

  // Per-phase timing summary (only when CCAKE_PROFILE is set; otherwise silent).
  if ( ccake::PhaseProfiler::instance().enabled() )
  {
    const double wall = sw.printTime();
    const double sps  = wall > 0.0 ? static_cast<double>(n_steps_ran)/wall : 0.0;
    formatted_output::summarize("Timesteps: " + to_string(n_steps_ran)
                                + "   (" + to_string(sps) + " steps/s)");
    ccake::PhaseProfiler::instance().report(std::cout);
  }

  // Write all buffered data to disk (once, at end)
  if (settingsPtr->particlization_enabled)
    outPtr->flush_freeze_out();
  if (settingsPtr->print_conservation_status)
    outPtr->flush_conservation(out_dir + "/conservation.dat");
  if (settingsPtr->check_causality)
    outPtr->flush_causality(out_dir + "/causality_summary.dat");
  // Final state of the BBMG partons, input for the fragmentation post-processor
  if (settingsPtr->jets_type != "disabled")
  {
    wsPtr->bbmg.copy_device_to_host_BBMG();
    outPtr->print_jet_final_to_txt(out_dir + "/jet_final.dat");
  }
}
