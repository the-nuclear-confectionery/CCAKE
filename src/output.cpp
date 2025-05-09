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
/// smoothed entropy density, extensive entropy, sigma_lab, sph_mass of entropy,
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

  vector<string> dataset_names = {"x", "y", "eta", "ux", "uy", "ueta",
                                  "T", "muB", "muS", "muQ",
																	"e", "s", "B", "S", "Q"};
  vector<string> dataset_units = {"fm", "fm", "dimensionless",
																	"dimensionless", "dimensionless", "1/fm",
																	"MeV", "MeV", "MeV", "MeV",
                                  "MeV/fm^3", "1/fm^3", "1/fm^3", "1/fm^3",
                                  "1/fm^3"};

  std::map<string,int> eos_map = {{"table",              0},
                                  {"tanh_conformal",     1},
                                  {"conformal",          2},
                                  {"conformal_diagonal", 3}};

	// cout << "dataset_names.size() = " << dataset_names.size() << endl;
	// cout << "dataset_units.size() = " << dataset_units.size() << endl;

  vector<vector<double> > data( dataset_names.size(),
                                vector<double>(systemPtr->particles.size()) );
  vector<int> eos_tags(systemPtr->particles.size());
  for (auto & p : systemPtr->particles)
  {
		switch (D)
		{
      case 1:
		    data[0][p.ID] = 0.0;
		    data[1][p.ID] = 0.0;
		    data[2][p.ID] = p.r(0);
		    data[3][p.ID] = 0.0;
		    data[4][p.ID] = 0.0;
		    data[5][p.ID] = p.hydro.u(0);
				break;
      case 2:
		    data[0][p.ID] = p.r(0);
		    data[1][p.ID] = p.r(1);
		    data[2][p.ID] = 0.0;
		    data[3][p.ID] = p.hydro.u(0);
		    data[4][p.ID] = p.hydro.u(1);
		    data[5][p.ID] = 0.0;
				break;
      case 3:
		    data[0][p.ID] = p.r(0);
		    data[1][p.ID] = p.r(1);
		    data[2][p.ID] = p.r(2);
		    data[3][p.ID] = p.hydro.u(0);
		    data[4][p.ID] = p.hydro.u(1);
		    data[5][p.ID] = p.hydro.u(2);
				break;
		}
    data[6][p.ID]  = p.T()*hbarc_MeVfm;
    data[7][p.ID]  = p.muB()*hbarc_MeVfm;
    data[8][p.ID]  = p.muS()*hbarc_MeVfm;
    data[9][p.ID]  = p.muQ()*hbarc_MeVfm;
    data[10][p.ID] = p.e()*hbarc_MeVfm;
    data[11][p.ID] = p.s();
    data[12][p.ID] = p.rhoB();
    data[13][p.ID] = p.rhoS();
    data[14][p.ID] = p.rhoQ();
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
       << result_cs2fzfluc(i) <<
       endl;
    count++;
  }
  formatted_output::detail("Printed " + std::to_string(count)
														+ " freeze-out particles.");
  formatted_output::detail("Current running total of "
														+ std::to_string(systemPtr->number_part_fo)
														+ " === " +  + "% freeze-out particles.");
	double total_pc = 100.0 * static_cast<double>(systemPtr->number_part_fo)
										/ static_cast<double>(systemPtr->number_part_fo_at_t0);
  formatted_output::detail("Currently " + std::to_string(total_pc)
														+ "% of all particles are frozen out.");
	int count_frozen_out_since_t0 = systemPtr->number_part_fo
																	- systemPtr->number_part_fo_at_t0;
	int total_not_frozen_out_at_t0 = systemPtr->n_particles
																		- systemPtr->number_part_fo_at_t0;
	double relative_pc = 100.0 * static_cast<double>(count_frozen_out_since_t0);
												/ static_cast<double>(total_not_frozen_out_at_t0);
  formatted_output::detail("Currently " + std::to_string(relative_pc)
														+ "% of all particles initially unfrozen out now "
															"frozen out.");


  formatted_output::detail("Currently " + std::to_string(total_pc)
														+ "% of all particles are frozen out.");


  FO.close();

  return;
}


/// @brief Prints the freeze-out particles to a text file.
/// @tparam D The dimensionality of the simulation.
/// @param[in] freeze_out Pointer to the freeze-out object.
/// @todo We need to understand the meaning of the quantities printed here.
template<>
void Output<2>::print_freeze_out(std::shared_ptr<FreezeOut<2>> freeze_out)
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
