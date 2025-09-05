#include "output.h"

using namespace constants;
using namespace ccake;


//Template instantiations
template class Output<1, EoM_default>;
template class Output<2, EoM_default>;
template class Output<3, EoM_default>;
template class Output<1, EoM_cartesian>;
template class Output<2, EoM_cartesian>;
template class Output<3, EoM_cartesian>;

// Constructors and destructors.

/// @brief Constructor for the Output class.
/// @details This constructor initializes the Output class with the pointers to
/// the Settings and SystemState objects. It also initializes the hdf5 file
/// if the hdf_evolution flag is set to true in the settings object.
/// @param[in] settingsPtr_in Pointer to the settings object.
/// @param[in] sys_in Pointer to the system state object.
/// @param[in] ws_in Pointer to the sph workstation object for bbmg.
template<unsigned int D, template<unsigned int> class TEOM>
Output<D,TEOM>::Output(std::shared_ptr<Settings> settingsPtr_in,
       std::shared_ptr<SystemState<D>> sys_in,
        std::shared_ptr<SPHWorkstation<D,TEOM>> ws_in)
    : settingsPtr(settingsPtr_in),
      systemPtr(sys_in),
      wsPtr(ws_in)
{
    std::vector<double> global_parameters_to_HDF = {
        settingsPtr->hT,
        settingsPtr->e_cutoff
    };
    std::vector<std::string> global_parameter_names_to_HDF = {
        "h",
        "e_cutoff"
    };
    output_directory = settingsPtr->results_directory;
    if (settingsPtr->txt_evolution)
        formatted_output::detail("WARNING: Output to txt enabled. Large files will be produced.");
    if (settingsPtr->hdf_evolution)
        hdf5_file.initialize(output_directory + "/system_state.h5",
                             global_parameters_to_HDF,
                             global_parameter_names_to_HDF);
}

template<unsigned int D, template<unsigned int> class TEOM> Output<D,TEOM>::~Output(){}


/// @brief Sets the path to the results directory.
/// \param[in] path_to_results_directory Path to the results directory.
template<unsigned int D, template<unsigned int> class TEOM>
void Output<D,TEOM>::set_results_directory( string path_to_results_directory )
{
  output_directory = path_to_results_directory;
}

//------------------------------------------------------------------------------
template<unsigned int D, template<unsigned int> class TEOM>
void Output<D,TEOM>::print_system_state()
{
  systemPtr->copy_device_to_host();
  //---------------------------------
  if (settingsPtr->txt_evolution)
    print_system_state_to_txt();

  //---------------------------------
  if (settingsPtr->hdf_evolution)
    print_system_state_to_HDF();

  if(settingsPtr->jet_evolution){
    if (systemPtr->n_particles == systemPtr->number_part_fo)
      {
      cout << systemPtr->n_particles << " is number of total particles; " << systemPtr->number_part_fo << " is number of part. frozen out. " << endl << "------------------------------------------------------" << endl;
    //cout << "printing random things from bbmg class" << bbmgPtr->jetInfo_host[1].T << endl;
      cout << "Copy device to host for jets executing..." << endl;
      wsPtr->bbmg.copy_device_to_host_BBMG();
      cout << "Print function for jets executing..." << endl;
      print_jet_state_to_txt();
      }
  }

  //---------------------------------
  // increment timestep index
  n_timesteps_output++;

  return;
}

/// @brief Print the system state to a text file.
/// @details This function prints the system state to a text file. There is a 
/// single initial line with the current time of the simulation. Then, each
/// particle is printed in a separate line. The columns are:
/// Particle ID, time, position, pressure, temperature, 
/// chemical potentials, energy density, baryon density, 
/// strangeness density, electric charge density, entropy density,
/// smoothed entropy density, extensive entropy, sigma_lab, sph_mass of entropy,
/// shear relaxation time, bigtheta, \f$\sqrt{\pi{\mu\nu}\pi^{\mu\nu}}\f$,
/// \f$\tau_{\pi}\Theta/\pi\f$, \f$\pi^{\tau\tau}\f$, \f$\pi^{xx}\f$,
/// \f$\pi^{yy}\f$, \f$\tau^2\pi^{\eta\eta}\f$, the fluid velocity, 
/// the Lorentz factor \f$\gamma\f$, the freeze-out flag, and the EOS name.
/// When relevant, quantities are converted to MeV and MeV/fm\f$^{3}\f$.
/// @tparam D The dimensionality of the simulation.
/// @todo: This needs to be updated to deal with other dimensions than 2+1D 
/// simulations, specially for the shear tensor.
template<unsigned int D, template<unsigned int> class TEOM>
void Output<D,TEOM>::print_system_state_to_txt()
{
  //if (n_timesteps_output > 30) {
  //  std::cout << "Terminating after the 30th timestep" << std::endl;
  //  exit(1);
  //}
  string outputfilename = output_directory + "/system_state_"
                          + std::to_string(n_timesteps_output) + ".dat";
  ofstream out( outputfilename.c_str() );
  out << systemPtr->t << "\n";
  int iParticle = 0;

  for ( auto & p : systemPtr->particles ){
    out << iParticle++ << " " //0
        << systemPtr->t << " " //1
//        << std::setw(8)
        << std::setprecision(6) << std::scientific;
        //2,3,4
        for (int idir=0; idir<D; idir++)
          out << p.r(idir) << " " ;
        for (int idir=D; idir<3; idir++)
          out << 0.0 << " ";
        out << p.p()*hbarc_GeVfm << " "       //5
        << p.T()*hbarc_GeVfm << " " //6
        << p.muB()*hbarc_GeVfm << " " //7
        << p.muS()*hbarc_GeVfm << " " //8
        << p.muQ()*hbarc_GeVfm << " " //9
        << p.e()*hbarc_GeVfm << " " //10
        << p.rhoB() << " " //11
        << p.rhoS() << " " //12
        << p.rhoQ() << " " //13
        << p.s() << " " //14
        << p.hydro.eta_pi << " " //15
        << p.hydro.zeta_Pi << " " //16
        << p.hydro.tau_Pi << " " //17
        << p.hydro.tau_pi << " " //18
        << p.hydro.theta << " "  //19
        << p.hydro.inverse_reynolds_shear << " " //20
        << p.hydro.inverse_reynolds_bulk << " " //21
        << p.hydro.shear_knudsen << " " //22
        << p.hydro.bulk_knudsen << " " //23
        << p.hydro.shv(0,0) << " " //24
        << p.hydro.shv(1,1) << " " //25
        << p.hydro.shv(2,2) << " " //26
        << p.hydro.shv(0,1) << " " //27
        << p.hydro.shv(0,2) << " " //28
        << p.hydro.shv(1,2) << " " //29
        << p.hydro.shv(1,3) << " " //30
        << p.hydro.shv(2,3) << " " //31
        << p.hydro.shv(3,3) << " "; //32
        //33,34, 35
        for (int idir=0; idir<D; idir++)
          out << p.hydro.u(idir) << " ";
        for (int idir=D; idir<3; idir++)
          out << 0.0 << " ";
        out << p.hydro.gamma << " " //36
        << p.Freeze << " " //37
        << p.get_current_eos_name() << " " //38
        << p.hydro.causality << " " << "\n";
  }

  out << std::flush;

  out.close();

  return;
}

/// @brief Outputs data from the uncoupled jet system to a dat file
/// @details This function outputs all necessary jet information for post processing. Each jet
/// is its own row, while each column is its own variable for that jet. The columns are as
/// follows: Index, time, BBMG line integral, initial density of each jet, angle the jet 
/// started, x position, and y position. All quantities are left in terms of fm and will be 
/// converted to MeV when necessary in post processing.
template<unsigned int D, template<unsigned int> class TEOM>
void Output<D,TEOM>::print_jet_state_to_txt()
{
  string jet_output_filename = output_directory + "/jet_state_" + std::to_string(n_timesteps_output) + ".dat";
  ofstream out( jet_output_filename.c_str() );
  // out << systemPtr->t << "\n";

  //cout << "(bbmgPtr->jetInfo).size() = " << (bbmgPtr->jetInfo).size() << endl;
  int iJet = 0;
  //cout << "(bbmgPtr->jetFreezeOut).size() = " << (bbmgPtr->jetFreezeOut).size() << endl;
  for ( auto & jets : wsPtr->bbmg.jetFreezeOut )//I think i wanna call the objects here jets, like how we have particles in the other function
       {
        out << iJet++ << " ";
        out << systemPtr->t << " ";
        out << jets.line_int << " ";
        out << jets.rho0 << " ";
        //out << jets.rho << " ";
        out << jets.PID << " ";
        out << jets.r[0] << " ";
        out << jets.r[1] << " " << "\n";
       }
  out << std::flush;
  out.close();
  return;
}

/// @brief Outputs data of the system state to an HDF5 file.
/// @tparam D The dimensionality of the simulation.
/// @details This function outputs the data of the system state to an HDF5 file.
/// The data is stored in datasets with the following names:
/// - x: x-coordinate of the particle.
/// - y: y-coordinate of the particle.
/// - T: Temperature of the particle.
/// - muB: Baryon chemical potential of the particle.
/// - muS: Strangeness chemical potential of the particle.
/// - muQ: Electric charge chemical potential of the particle.
/// - e: Energy density of the particle.
/// - s: Entropy density of the particle.
/// - B: Baryon density of the particle.
/// - S: Strangeness density of the particle.
/// - Q: Electric charge density of the particle.
/// The data is stored in the datasets in the order of the particles.
/// The units of the data are stored in the dataset attributes.
/// @todo: This needs to be adapated for other dimensions than 2+1D simulations.
template<unsigned int D, template<unsigned int> class TEOM>
void Output<D,TEOM>::print_system_state_to_HDF()
{
  // get width from maximum possible number of timesteps
  const int width = ceil(log10(ceil(settingsPtr->max_tau/settingsPtr->dt)));

  vector<string> dataset_names = {"x", "y", "T", "muB", "muS", "muQ",
                                  "e", "s", "B", "S", "Q"};
  vector<string> dataset_units = {"fm", "fm", "MeV", "MeV", "MeV", "MeV",
                                  "MeV/fm^3", "1/fm^3", "1/fm^3", "1/fm^3",
                                  "1/fm^3"};

  std::map<string,int> eos_map = {{"table",              0},
                                  {"tanh_conformal",     1},
                                  {"conformal",          2},
                                  {"conformal_diagonal", 3}};

  vector<vector<double> > data( dataset_names.size(),
                                vector<double>(systemPtr->particles.size()) );
  vector<int> eos_tags(systemPtr->particles.size());
  for (auto & p : systemPtr->particles)
  {
    data[0][p.ID]  = p.r(0);
    data[1][p.ID]  = p.r(1);
    data[2][p.ID]  = p.T()*hbarc_MeVfm;
    data[3][p.ID]  = p.muB()*hbarc_MeVfm;
    data[4][p.ID]  = p.muS()*hbarc_MeVfm;
    data[5][p.ID]  = p.muQ()*hbarc_MeVfm;
    data[6][p.ID]  = p.e()*hbarc_MeVfm;
    data[7][p.ID]  = p.s();
    data[8][p.ID]  = p.rhoB();
    data[9][p.ID]  = p.rhoS();
    data[10][p.ID] = p.rhoQ();
    eos_tags[p.ID] = eos_map[ p.get_current_eos_name() ];
  }

  vector<string> parameter_names = { "Time", "e_2_X", "e_2_P" };
  //vector<double> parameters      = { systemPtr->t, systemPtr->e_2_X.back(),
  //                                                 systemPtr->e_2_P.back() };
  vector<double> parameters      = { systemPtr->t, -10, -10 }; //Disabled for now

  hdf5_file.output_dataset( dataset_names, dataset_units, data, width,
                            n_timesteps_output, eos_tags,
                            parameters, parameter_names );

  return;
}

/// @brief Prints the conservation status of the system.
/// @tparam D The dimensionality of the simulation.
template<unsigned int D, template<unsigned int> class TEOM>
void Output<D,TEOM>::print_conservation_status()
{   
    double Bdiff = (systemPtr->Btotal - systemPtr->Btotal0)/ systemPtr->Btotal0 * 100.0;
    double Sdiff = (systemPtr->Stotal - systemPtr->Stotal0)/ systemPtr->Stotal0 * 100.0;
    double Qdiff = (systemPtr->Qtotal - systemPtr->Qtotal0)/ systemPtr->Qtotal0 * 100.0;
    stringstream ss;
    //ss  << "t = "
    //    << systemPtr->t      << ": " << scientific        << setw(10)
    //    << systemPtr->Eloss  << " "  << systemPtr->S      << " "
    //    << systemPtr->Btotal << " "  << systemPtr->Stotal << " "
    //    << systemPtr->Qtotal << defaultfloat;
    ss  << "t = "
        << systemPtr->t      << ": " << std::scientific << std::setprecision(6)
        << systemPtr->Eloss  << " "  << systemPtr->S      << " "
        << Bdiff  << " "  << Sdiff << " "
        << Qdiff  << " "  << std::defaultfloat;
    formatted_output::summarize(ss.str());
    stringstream sss;
    sss << "Btotal = "
        << systemPtr->Btotal << " "
        << "Stotal = "
        << systemPtr->Stotal << " "
        << "Qtotal = "
        << systemPtr->Qtotal << " " << std::defaultfloat;
    formatted_output::summarize(sss.str());

}


/// @brief Prints the freeze-out particles to a text file.
/// @tparam D The dimensionality of the simulation.
/// @param[in] freeze_out Pointer to the freeze-out object.
/// @todo We need to understand the meaning of the quantities printed here.
template<unsigned int D, template<unsigned int> class TEOM>
void Output<D,TEOM>::print_freeze_out(std::shared_ptr<FreezeOut<D>> freeze_out, double B, double Q, double S)
{
  string outputfilename = output_directory + "/freeze_out.dat";
  ofstream FO( outputfilename.c_str(), ios::app );


  //bool write_header = false;
  //struct stat file_stat;
  //if (stat(outputfilename.c_str(), &file_stat) != 0 || file_stat.st_size == 0) {
  //  write_header = true;
  //}
//
  ////print total B,S,Q of the system
  //if (write_header) {
  //  FO << "# " << B << " "
  //             << S << " "
  //             << Q << "\n";
  //}

  //Copy data to host
  auto FOResults = Cabana::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                                        freeze_out->results);
  int count=0;
  FRZ_RESULTS_VIEW(result_, FOResults)
  for (int i = 0; i < FOResults.size(); i++){
    if (!result_print(i)) continue;
    FO << result_divEener(i) << " ";
    for(int idir = 0; idir < D; idir++)
      FO << result_divE(i, idir) << " ";
    FO << result_gsub(i) << " ";
    for(int idir = 0; idir < D; idir++)
      FO << result_uout(i, idir) << " ";
    FO << result_swsub(i) << " "
       << result_bulksub(i) << " "
       << result_shearsub(i, 0,0) << " "
       << result_shearsub(i, 1,1) << " "
       << result_shearsub(i, 2,2) << " "
       << result_shear33sub(i) << " "
       << result_shearsub(i,1,2) << " "
       << result_shearsub(i,1,3) << " "
       << result_shearsub(i,2,3) << " "
       << result_tlist(i) << " ";
    for(int idir = 0; idir < D; idir++)
        FO << result_rsub(i, idir) << " ";
    FO << result_sFO(i) << " "
       << result_Efluc(i) << " "
       << result_Tfluc(i) << " "
       << result_muBfluc(i) << " "
       << result_muSfluc(i) << " "
       << result_muQfluc(i) << " "
       << result_wfzfluc(i) << " "
       << result_cs2fzfluc(i) <<// " "
       //<< result_rhoBfluc(i) << " "
       //<< result_rhoSfluc(i) << " "
       //<< result_rhoQfluc(i) << 
       endl;
    count++;
  }
  formatted_output::detail("Printed " + std::to_string(count) + " freeze-out particles.");
  FO.close();

  return;
}


/// @brief Prints the freeze-out particles to a text file.
/// @tparam D The dimensionality of the simulation.
/// @param[in] freeze_out Pointer to the freeze-out object.
/// @todo We need to understand the meaning of the quantities printed here.
template<>
void Output<2,EoM_default>::print_freeze_out(std::shared_ptr<FreezeOut<2>> freeze_out, double B, double Q, double S)
{
  string outputfilename = output_directory + "/freeze_out.dat";
  ofstream FO( outputfilename.c_str(), ios::app );

  //Copy data to host
  auto FOResults = Cabana::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                                        freeze_out->results);

  bool write_header = false;
  struct stat file_stat;
  if (stat(outputfilename.c_str(), &file_stat) != 0 || file_stat.st_size == 0) {
    write_header = true;
  }

  //print total B,S,Q of the system
  //if (write_header) {
  //  FO << "# " << B << " "
  //             << S << " "
  //             << Q << "\n";
  //}

  int count=0;
  FRZ_RESULTS_VIEW(result_, FOResults)
  for (int i = 0; i < FOResults.size(); i++){
    if (!result_print(i)) continue;
    FO << result_divEener(i) << " ";
    for(int idir = 0; idir < 2; idir++)
      FO << result_divE(i, idir) << " ";
    FO << result_gsub(i) << " ";
    for(int idir = 0; idir < 2; idir++)
      FO << result_uout(i, idir) << " ";
    FO << result_swsub(i) << " "
       << result_bulksub(i) << " "
       << result_shearsub(i, 0,0) << " "
       << result_shearsub(i, 1,1) << " "
       << result_shearsub(i, 2,2) << " "
       << result_shear33sub(i) << " "
       << result_shearsub(i,1,2) << " "
       << result_tlist(i) << " ";
    for(int idir = 0; idir < 2; idir++)
        FO << result_rsub(i, idir) << " ";
    FO << result_sFO(i) << " "
       << result_Efluc(i) << " "
       << result_Tfluc(i) << " "
       << result_muBfluc(i) << " "
       << result_muSfluc(i) << " "
       << result_muQfluc(i) << " "
       << result_wfzfluc(i) << " "
       << result_cs2fzfluc(i) <<// " "
        //<< result_rhoBfluc(i) << " "
        //<< result_rhoSfluc(i) << " "
        //<< result_rhoQfluc(i) <<
       endl;
    count++;
  }
  formatted_output::detail("Printed " + std::to_string(count) + " freeze-out particles.");
  FO.close();

  return;
}

template<>
void Output<2,EoM_cartesian>::print_freeze_out(std::shared_ptr<FreezeOut<2>> freeze_out, double B, double S, double Q)
{
  string outputfilename = output_directory + "/freeze_out.dat";
  ofstream FO( outputfilename.c_str(), ios::app );

  //Copy data to host
  auto FOResults = Cabana::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                                        freeze_out->results);
  int count=0;
  FRZ_RESULTS_VIEW(result_, FOResults)
  for (int i = 0; i < FOResults.size(); i++){
    if (!result_print(i)) continue;
    FO << result_divEener(i) << " ";
    for(int idir = 0; idir < 2; idir++)
      FO << result_divE(i, idir) << " ";
    FO << result_gsub(i) << " ";
    for(int idir = 0; idir < 2; idir++)
      FO << result_uout(i, idir) << " ";
    FO << result_swsub(i) << " "
       << result_bulksub(i) << " "
       << result_shearsub(i, 0,0) << " "
       << result_shearsub(i, 1,1) << " "
       << result_shearsub(i, 2,2) << " "
       << result_shear33sub(i) << " "
       << result_shearsub(i,1,2) << " "
       << result_tlist(i) << " ";
    for(int idir = 0; idir < 2; idir++)
        FO << result_rsub(i, idir) << " ";
    FO << result_sFO(i) << " "
       << result_Efluc(i) << " "
       << result_Tfluc(i) << " "
       << result_muBfluc(i) << " "
       << result_muSfluc(i) << " "
       << result_muQfluc(i) << " "
       << result_wfzfluc(i) << " "
       << result_cs2fzfluc(i) <<
       endl;
    count++;
  }
  formatted_output::detail("Printed " + std::to_string(count) + " freeze-out particles.");
  FO.close();

  return;
}
