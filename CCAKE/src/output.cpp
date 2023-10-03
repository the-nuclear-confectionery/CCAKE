#include "output.h"

using namespace constants;
using namespace ccake;


//Template instantiations
template class Output<1>;
template class Output<2>;
template class Output<3>;

// Constructors and destructors.

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
  hdf5_file.initialize( output_directory + "/system_state.h5",
                        global_parameters_to_HDF,
                        global_parameter_names_to_HDF );
}

template<unsigned int D> Output<D>::~Output(){}


/// \brief Sets the path to the results directory.
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
  if (settingsPtr->printing_to_txt)
    print_system_state_to_txt();

  //---------------------------------
  if (settingsPtr->printing_to_HDF)
    print_system_state_to_HDF();

  //---------------------------------
  //print_freeze_out();

  //---------------------------------
  // increment timestep index
  n_timesteps_output++;

  return;
}

//------------------------------------------------------------------------------
template<unsigned int D>
void Output<D>::print_system_state_to_txt()
{
  string outputfilename = output_directory + "/system_state_"
                          + std::to_string(n_timesteps_output) + ".dat";
  ofstream out( outputfilename.c_str() );

  out << systemPtr->t << "\n";
  int iParticle = 0;
  
  for ( auto & p : systemPtr->particles )
    out << iParticle++ << " "
        << systemPtr->t << " "
        << p.r << " "
        << p.p() << " "
        << p.T()*hbarc_MeVfm << " "
        << p.muB()*hbarc_MeVfm << " "
        << p.muS()*hbarc_MeVfm << " "
        << p.muQ()*hbarc_MeVfm << " "
        << p.e()*hbarc_MeVfm << "       " //10
        << p.rhoB() << " "
        << p.rhoS() << " "
        << p.rhoQ() << " "
        << p.s() << " "
        << p.smoothed.s/(p.hydro.gamma*systemPtr->t) << " "
        << p.specific.s << " "
        << p.hydro.sigma << " "
        << p.norm_spec.s << " "
        << p.hydro.stauRelax << " "
        << p.hydro.bigtheta << "       "  //20
        << sqrt(     p.hydro.shv(0,0)*p.hydro.shv(0,0)
                -2.0*p.hydro.shv(0,1)*p.hydro.shv(0,1)
                -2.0*p.hydro.shv(0,2)*p.hydro.shv(0,2)
                +    p.hydro.shv(1,1)*p.hydro.shv(1,1)
                +    p.hydro.shv(2,2)*p.hydro.shv(2,2)
                +2.0*p.hydro.shv(1,2)*p.hydro.shv(1,2)
                +pow(systemPtr->t,4.0)*p.hydro.shv33*p.hydro.shv33 ) << " "
        << p.hydro.stauRelax/systemPtr->t * p.hydro.bigtheta << " "
        << p.hydro.shv(0,0) << " "
        << p.hydro.shv(1,1) << " "
        << p.hydro.shv(2,2) << " "
        << p.hydro.shv(1,2) << " "
        << pow(systemPtr->t,2.0)*p.hydro.shv33 << " "
        << p.hydro.u(0)/p.hydro.gamma << " "  //28
        << p.hydro.u(1)/p.hydro.gamma << " "
        << p.hydro.gamma << "       "
        << p.Freeze << " "
        << p.get_current_eos_name() << "\n";

  out << std::flush;

  out.close();

  return;
}

//------------------------------------------------------------------------------
template<unsigned int D>
void Output<D>::print_system_state_to_HDF()
{
  // get width from maximum possible number of timesteps
  const int width = ceil(log10(ceil(settingsPtr->tend/settingsPtr->dt)));

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


//==============================================================================
/*template<unsigned int D>
void Output<D>::print_freeze_out()
{
  string outputfilename = output_directory + "/freeze_out_"
                          + std::to_string(n_timesteps_output) + ".dat";
  ofstream FO( outputfilename.c_str(), ios::out | ios::app );

  auto & freeze_out = wsPtr->freeze_out;

  if ( freeze_out.divTtemp.size() != freeze_out.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( freeze_out.divT.size() != freeze_out.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( freeze_out.gsub.size() != freeze_out.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( freeze_out.uout.size() != freeze_out.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( freeze_out.swsub.size() != freeze_out.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( freeze_out.bulksub.size() != freeze_out.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( freeze_out.shearsub.size() != freeze_out.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( freeze_out.shear33sub.size() != freeze_out.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( freeze_out.tlist.size() != freeze_out.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( freeze_out.rsub.size() != freeze_out.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( freeze_out.sFO.size() != freeze_out.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( freeze_out.Tfluc.size() != freeze_out.cf ) cout << "WARNING: WRONG SIZE" << endl;

  for (int i = 0; i < freeze_out.cf; i++)
    FO << freeze_out.divTtemp[i] << " " << freeze_out.divT[i] << " "
        << freeze_out.gsub[i] << " " << freeze_out.uout[i] << " "
        << freeze_out.swsub[i] << " " << freeze_out.bulksub[i] << " "
        << freeze_out.shearsub[i](0,0) << " "
        << freeze_out.shearsub[i](1,1) << " "
        << freeze_out.shearsub[i](2,2) << " "
        << freeze_out.shear33sub[i] << " "
        << freeze_out.shearsub[i](1,2) << " "
        << freeze_out.tlist[i] << " " << freeze_out.rsub[i] << " "
        << freeze_out.sFO[i] << " " << freeze_out.Tfluc[i] << endl;

  FO.close();

  return;
}*/
