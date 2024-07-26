#include "output.h"

using namespace constants;
using namespace ccake;


//Template instantiations
template class Output<1>;
template class Output<2>;
template class Output<3>;

// Constructors and destructors.

/// @brief Constructor for the Output class.
/// @details This constructor initializes the Output class with the pointers to
/// the Settings and SystemState objects. It also initializes the hdf5 file
/// if the hdf_evolution flag is set to true in the settings object.
/// @param[in] settingsPtr_in Pointer to the settings object.
/// @param[in] sys_in Pointer to the system state object.
template<unsigned int D>
Output<D>::Output( std::shared_ptr<Settings> settingsPtr_in,
                   std::shared_ptr<SystemState<D>> sys_in ):
                   settingsPtr(settingsPtr_in), systemPtr(sys_in)
{
  vector<double> global_parameters_to_HDF
                  = vector<double>({ settingsPtr->hT,
                                     settingsPtr->e_cutoff });
  vector<string> global_parameter_names_to_HDF
                  = vector<string>({ "h",
                                     "e_cutoff" });
  output_directory = settingsPtr->results_directory;
  if (settingsPtr->txt_evolution)
    formatted_output::detail("WARNING: Output to txt enabled. Large files will be produced.");
  if (settingsPtr->hdf_evolution)
    hdf5_file.initialize( output_directory + "/system_state.h5",
                        global_parameters_to_HDF,
                        global_parameter_names_to_HDF );
}

template<unsigned int D> Output<D>::~Output(){}


/// @brief Sets the path to the results directory.
/// \param[in] path_to_results_directory Path to the results directory.
template<unsigned int D>
void Output<D>::set_results_directory( string path_to_results_directory )
{
  output_directory = path_to_results_directory;
}

//------------------------------------------------------------------------------
template<unsigned int D>
void Output<D>::print_system_state()
{
  systemPtr->copy_device_to_host();
  //---------------------------------
  if (settingsPtr->txt_evolution)
    print_system_state_to_txt();

  //---------------------------------
  if (settingsPtr->hdf_evolution)
    print_system_state_to_HDF();

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
/// smoothed entropy density, specific entropy, sigma, norm_spec of entropy,
/// shear relaxation time, bigtheta, \f$\sqrt{\pi{\mu\nu}\pi^{\mu\nu}}\f$,
/// \f$\tau_{\pi}\Theta/\pi\f$, \f$\pi^{\tau\tau}\f$, \f$\pi^{xx}\f$,
/// \f$\pi^{yy}\f$, \f$\tau^2\pi^{\eta\eta}\f$, the fluid velocity, 
/// the Lorentz factor \f$\gamma\f$, the freeze-out flag, and the EOS name.
/// When relevant, quantities are converted to MeV and MeV/fm\f$^{3}\f$.
/// @tparam D The dimensionality of the simulation.
/// @todo: This needs to be updated to deal with other dimensions than 2+1D 
/// simulations, specially for the shear tensor.
template<unsigned int D>
void Output<D>::print_system_state_to_txt()
{
  if (n_timesteps_output > 30) {
    std::cout << "Terminating after the 30th timestep" << std::endl;
    exit(1);
  }
  string outputfilename = output_directory + "/system_state_"
                          + std::to_string(n_timesteps_output) + ".dat";
  ofstream out( outputfilename.c_str() );
  out << systemPtr->t << "\n";
  int iParticle = 0;

  for ( auto & p : systemPtr->particles ){
    out << iParticle++ << " " //0
        << systemPtr->t << " " //1
//        << std::setw(8)
        << std::setprecision(6) << std::scientific
        << p.r          //2,3 (in 2D)
        << p.p() << " "       //4
        << p.T()*hbarc_MeVfm << " " //5
        << p.muB()*hbarc_MeVfm << " " //6
        << p.muS()*hbarc_MeVfm << " " //7
        << p.muQ()*hbarc_MeVfm << " " //8
        << p.e()*hbarc_MeVfm << " " //9
        << p.rhoB() << " " //10
        << p.rhoS() << " " //11
        << p.rhoQ() << " " //12
        << p.s() << " " //13
        << p.smoothed.s/(p.hydro.gamma*systemPtr->t) << " " //14
        << p.specific.s << " " //15
        << p.hydro.sigma << " " //16
        << p.norm_spec.s << " " //17
        << p.hydro.stauRelax << " " //18
        << p.hydro.bigtheta << " "  //19
        << sqrt(     p.hydro.shv(0,0)*p.hydro.shv(0,0)
                -2.0*p.hydro.shv(0,1)*p.hydro.shv(0,1)
                -2.0*p.hydro.shv(0,2)*p.hydro.shv(0,2)
                +    p.hydro.shv(1,1)*p.hydro.shv(1,1)
                +    p.hydro.shv(2,2)*p.hydro.shv(2,2)
                +2.0*p.hydro.shv(1,2)*p.hydro.shv(1,2)
                +pow(systemPtr->t,4.0)*p.hydro.shv33*p.hydro.shv33 ) << " " //20
        << p.hydro.stauRelax/systemPtr->t * p.hydro.bigtheta << " " //21
        << p.hydro.shv(0,0) << " " //22
        << p.hydro.shv(1,1) << " " //23
        << p.hydro.shv(2,2) << " " //24
        << p.hydro.shv(1,2) << " " //25
        << pow(systemPtr->t,2.0)*p.hydro.shv33 << " "; //26
        for (int idir=0; idir<D; idir++)
          out << p.hydro.u(idir)/p.hydro.gamma << " "; //27,28 (in 2+1)
        out << p.hydro.gamma << " " //29
        << p.Freeze << " " //30
        << p.get_current_eos_name() << "\n"; //31
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
template<unsigned int D>
void Output<D>::print_system_state_to_HDF()
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
template<unsigned int D>
void Output<D>::print_conservation_status()
{   
    stringstream ss;
    ss  << "t = "
        << systemPtr->t      << ": " << scientific        << setw(10)
        << systemPtr->Eloss  << " "  << systemPtr->S      << " "
        << systemPtr->Btotal << " "  << systemPtr->Stotal << " "
        << systemPtr->Qtotal << defaultfloat;
    formatted_output::summarize(ss.str());
}


/// @brief Prints the freeze-out particles to a text file.
/// @tparam D The dimensionality of the simulation.
/// @param[in] freeze_out Pointer to the freeze-out object.
/// @todo We need to understand the meaning of the quantities printed here.
template<unsigned int D>
void Output<D>::print_freeze_out(std::shared_ptr<FreezeOut<D>> freeze_out)
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
       << result_cs2fzfluc(i) <<
       endl;
    count++;
  }
  formatted_output::detail("Printed " + std::to_string(count) + " freeze-out particles.");
  FO.close();

  return;
}
