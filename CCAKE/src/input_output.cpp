#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "../include/constants.h"
#include "../include/defaults.h"
#include "../include/formatted_output.h"
#include "../include/input_output.h"
#include "../include/particle.h"
#include "../include/mathdef.h"
#include "../include/sph_workstation.h"
#include "../include/vector.h"

using namespace constants;

using std::endl;
using std::flush;
using std::string;
using std::to_string;

// Constructors and destructors.
InputOutput::InputOutput(){}
InputOutput::~InputOutput(){}

//------------------------------------------------------------------------------
void InputOutput::set_EquationOfStatePtr( EquationOfState * eosPtr_in )
{
  eosPtr = eosPtr_in;
}

//------------------------------------------------------------------------------
void InputOutput::set_SettingsPtr( Settings * settingsPtr_in )
{
  settingsPtr = settingsPtr_in;
}

//------------------------------------------------------------------------------
void InputOutput::set_SystemStatePtr( SystemState * systemPtr_in )
{
  systemPtr = systemPtr_in;
}

//------------------------------------------------------------------------------
void InputOutput::set_SPHWorkstationPtr( SPHWorkstation * wsPtr_in )
{
  wsPtr = wsPtr_in;
}

//------------------------------------------------------------------------------
void InputOutput::set_results_directory( string path_to_results_directory )
{
  output_directory = path_to_results_directory;
}

//------------------------------------------------------------------------------
void InputOutput::load_settings_file( string path_to_settings_file )
{
  formatted_output::report("Reading in input parameter settings");

  string Param_file = path_to_settings_file;
  ifstream infile( Param_file.c_str() );
  if (infile.is_open())
  {
    string line;
    string field = "";
    string value = "";
    
    //--------------------------------------------------------------------------
    // set default values first
    setting_map values = parameter_settings::get_defaults();

    //--------------------------------------------------------------------------
    // update settings from input file
    while ( getline (infile, line) )
    {
      // skip whitespace lines
      if ( line.find_first_not_of(" \t\n\v\f\r") == std::string::npos )
        continue;

      istringstream iss(line);
      iss >> field >> value;

      // skip commented lines
      if ( field.front() == comment_character )
        continue;

      remove_char(field, ':');
      set_value(values, field, value);
    }

    //--------------------------------------------------------------------------
    // store final settings
    settingsPtr->IC_type                =      get_value(values, "ICtype");
    settingsPtr->IC_option              =      get_value(values, "ICoption");
    settingsPtr->IC_file                =      get_value(values, "ICfile");

    settingsPtr->h                      = stod(get_value(values, "h"));
    settingsPtr->dt                     = stod(get_value(values, "dt"));
    settingsPtr->t0                     = stod(get_value(values, "t0"));

    settingsPtr->e_cutoff               = stod(get_value(values, "e_cutoff"))/hbarc_GeVfm;

    settingsPtr->EoS_type               =      get_value(values, "EoS_Type");
    settingsPtr->EoS_option             =      get_value(values, "EoS_Option");

    settingsPtr->etaMode                =      get_value(values, "etaMode");
    settingsPtr->constant_eta_over_s    = stod(get_value(values, "constant_eta_over_s"));
    settingsPtr->shearRelaxMode         =      get_value(values, "shearRelaxMode");
    settingsPtr->zetaMode               =      get_value(values, "zetaMode");
    settingsPtr->constant_zeta_over_s   = stod(get_value(values, "constant_zeta_over_s"));
    settingsPtr->cs2_dependent_zeta_A   = stod(get_value(values, "cs2_dependent_zeta_A"));
    settingsPtr->cs2_dependent_zeta_p   = stod(get_value(values, "cs2_dependent_zeta_p"));
    settingsPtr->bulkRelaxMode          =      get_value(values, "bulkRelaxMode");

    settingsPtr->Freeze_Out_Temperature = stod(get_value(values, "freezeoutT"))/hbarc_MeVfm;
    settingsPtr->Freeze_Out_Type        =      get_value(values, "freezeout");


    //--------------------------------------------------------------------------
    // check set values
    formatted_output::update("Check parameter settings:");
    for ( auto & value : values )
      formatted_output::detail( "set " + value.first + " == " + value.second );

    infile.close();
  }
  else
  {
    std::cerr << "File " << Param_file << " could not be opened!\n";
    abort();
  }

  // make sure that chosen settings make sense
  settingsPtr->check_consistency();


  // set particles to print
  settingsPtr->particles_to_print = vector<int>({});


  // set up HDF5 output file here
  vector<double> global_parameters_to_HDF
                  = vector<double>({ settingsPtr->h,
                                     settingsPtr->e_cutoff });
  vector<string> global_parameter_names_to_HDF
                  = vector<string>({ "h",
                                     "e_cutoff" });
  hdf5_file.initialize( output_directory + "/system_state.h5",
                        global_parameters_to_HDF,
                        global_parameter_names_to_HDF );


  return;
}

//------------------------------------------------------------------------------
void InputOutput::set_EoS_type()
{
  string EoS_type           = settingsPtr->EoS_type;
  string EoS_option         = settingsPtr->EoS_option;
  string EoS_files_location = "EoS/" + EoS_type + "/" + EoS_option;
  string densities          = EoS_files_location + "/densities.dat";
  string derivatives        = EoS_files_location + "/derivatives.dat";

  if (EoS_option != "default")
  {
    std::cerr << "EoS option not recognized for " << EoS_type
              << ", now exiting.\n";
    abort();
  }

  eosPtr->quantity_file = densities;
  eosPtr->deriv_file    = derivatives;

  return;
}

//------------------------------------------------------------------------------
void InputOutput::read_in_initial_conditions()
{
  formatted_output::report("Reading in initial conditions for hydrodynamics");

  string initial_condition_type = settingsPtr->IC_type;
  formatted_output::update("Initial conditions type: " + settingsPtr->IC_type);

  int total_header_lines;
  string IC_file = settingsPtr->IC_file;

  //----------------------------------------------------------------------------
  if (initial_condition_type == "ICCING")
  {
    total_header_lines = 1;

    ifstream infile(IC_file.c_str());
    formatted_output::update("Initial conditions file: " + IC_file);
    if (infile.is_open())
    {
      string line;
      int count_header_lines = 0;
      int count_file_lines   = 0;
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

          Particle p;
          p.r(0)       = x;
          p.r(1)       = y;
          p.input.e    = e;
          p.input.rhoB = rhoB;
          p.input.rhoS = rhoS;
          p.input.rhoQ = rhoQ;
          p.hydro.u(0) = ux;
          p.hydro.u(1) = uy;

          systemPtr->particles.push_back( p );
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
  //----------------------------------------------------------------------------
  else if (initial_condition_type == "Freestream")
  {
    total_header_lines = 1;
    settingsPtr->initializing_with_full_Tmunu = true;

    ifstream infile(IC_file.c_str());
    formatted_output::update("Initial conditions file: " + IC_file);

    if (infile.is_open())
    {
      string line;
      int count_header_lines = 0;
      int count_file_lines   = 0;
      double x = 0.0, y = 0.0, e = 0.0, varsigma = 0.0, u1 = 0.0, u2 = 0.0;

      double pi00 = 0.0, pi01 = 0.0, pi02 = 0.0, pi10 = 0.0, pi11 = 0.0,
             pi12 = 0.0, pi20 = 0.0, pi21 = 0.0, pi22 = 0.0, pi33 = 0.0;
      double ignore = 0.0, stepX = 0.0, stepY = 0.0;

      while (getline (infile, line))
      {
        istringstream iss(line);
        if (count_header_lines < total_header_lines)
        {
          settingsPtr->headers.push_back(line);

          iss >> ignore >> stepX >> stepY;

          settingsPtr->stepx = stepX;
          settingsPtr->stepy = stepY;
          count_header_lines++;
        }
        else
        {
          iss >> x >> y >> e >> varsigma
              >> u1 >> u2
              >> pi00 >> pi01 >> pi02
              >> pi10 >> pi11 >> pi12
              >> pi20 >> pi21 >> pi22;

          e        /= hbarc_GeVfm;
          varsigma /= hbarc_GeVfm;
          pi00     /= hbarc_GeVfm;
          pi11     /= hbarc_GeVfm;
          pi12     /= hbarc_GeVfm;
          pi22     /= hbarc_GeVfm;

          pi33      = (pi00 - pi11 - pi22)/(settingsPtr->t0*settingsPtr->t0);

          // initialize particle object
          Particle p;
          p.r(0)           = x;
          p.r(1)           = y;
          p.input.e        = e;
          p.hydro.varsigma = varsigma;
          p.hydro.u(0)     = u1;
          p.hydro.u(1)     = u2;
          p.hydro.shv(0,0) = 0.0;
          p.hydro.shv(0,1) = 0.0;
          p.hydro.shv(0,2) = 0.0;
          p.hydro.shv(1,0) = 0.0;
          p.hydro.shv(1,1) = pi11;
          p.hydro.shv(1,2) = pi12;
          p.hydro.shv(2,0) = 0.0;
          p.hydro.shv(2,1) = pi12;
          p.hydro.shv(2,2) = pi22;
          p.hydro.shv33    = pi33;

          systemPtr->particles.push_back( p );
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
  //----------------------------------------------------------------------------
  else if (initial_condition_type == "Gubser")
  {
    // choose initial coordinate system
    const bool BSQmode = static_cast<bool>( settingsPtr->IC_option == "BSQ" );

    // initial time
    const double tau0 = settingsPtr->t0;

    // set Gubser profile parameters
    const double q     = 1.0;                   // 1/fm
    const double e0    = 1.0;                   // 1/fm^4
    const double rhoB0 = (BSQmode) ? 0.5 : 0.0; // 1/fm^3
    const double rhoQ0 = (BSQmode) ? 0.5 : 0.0; // 1/fm^3
    const double rhoS0 = (BSQmode) ? 0.5 : 0.0; // 1/fm^3

    // GRID GENERATION IN CARTESIAN COORDINATES
    // set grid step size for test
    const double TINY  = 1e-10;
    const double dx    = 0.025, dy = 0.025;
    settingsPtr->stepx = dx;
    settingsPtr->stepy = dy;
    const double xmin  = -5.0, xmax = 5.0+dx*TINY;
    const double ymin  = -5.0, ymax = 5.0+dy*TINY;

    // generate initial profile in (r,phi)
    double q2 = q*q, q4 = q2*q2, t2 = tau0*tau0, t3 = t2*tau0, t4 = t3*tau0;
    for ( double x = xmin; x <= xmax; x += dx )
    for ( double y = ymin; y <= ymax; y += dy )
    {
      double r         = sqrt(x*x+y*y);
      double r2        = r*r;
      double arg       = 1.0 + 2.0*q2*(t2+r2) + q4*(t2-r2)*(t2-r2);

      double eLocal    = (e0/t4)*pow(2.0*q*tau0, 8.0/3.0) / pow(arg, 4.0/3.0);
      double rhoBLocal = (rhoB0/t3)*4.0*q2*t2/arg;
      double rhoQLocal = (rhoQ0/t3)*4.0*q2*t2/arg;
      double rhoSLocal = (rhoS0/t3)*4.0*q2*t2/arg;

      double vr = 2.0*q2*tau0*r/(1+q2*t2+q2*r2);
      double gammar = 1.0/sqrt(1.0-vr*vr);

      double phi = atan2(y, x);
      double cphi = cos(phi), sphi = sin(phi);
      double ux = gammar*vr*cphi, uy = gammar*vr*sphi;

      Particle p;
      p.r(0)       = x;
      p.r(1)       = y;
      p.input.e    = eLocal;
      p.input.rhoB = rhoBLocal;
      p.input.rhoS = rhoSLocal;
      p.input.rhoQ = rhoQLocal;
      p.hydro.u(0) = ux;
      p.hydro.u(1) = uy;

      systemPtr->particles.push_back( p );
    }
    
  }
  //----------------------------------------------------------------------------
  else if (initial_condition_type == "Gubser_with_shear")
  { // NOTA BENE: THIS MODE MUST BE LOADED FROM FILE
    // choose initial coordinate system

    // initial time
    const double tau0 = settingsPtr->t0;

    // set Gubser profile parameters
    const double q     = 1.0; // 1/fm
    const double e0    = 9126.0*pi*pi/3125.0; // 1/fm^4
                              // use this normalization to compare with semi-
                              // analytic calculation in Phys. Rev. C 91, 014903
    const double rhoB0 = 0.0; // 1/fm^3
    const double rhoQ0 = 0.0; // 1/fm^3
    const double rhoS0 = 0.0; // 1/fm^3

    // GRID GENERATION IN CARTESIAN COORDINATES
    // set grid step size for test
    const double TINY  = 1e-10;
    const double dx    = 0.05, dy = 0.05;
    settingsPtr->stepx = dx;
    settingsPtr->stepy = dy;
    const double xmin  = -5.0, xmax = 5.0+dx*TINY;
    const double ymin  = -5.0, ymax = 5.0+dy*TINY;

    // sets the scale for conformal EoS
    const double Nc    = 3.0; // three colors
    const double Nf    = 2.5; // u+d massless, s 'half massless'
    const double cpLoc = pi*pi*(2.0*(Nc*Nc-1.0)+(7.0/2.0)*Nc*Nf)/90.0;

    // load input file
    string inputfilename = "./misc/Gubser_checks/ac/Initial_Profile_tau=1fm.dat";
    //cout << "Reading in Gubser initial profile from " << inputfilename << endl;
    formatted_output::update("Initial conditions file: " + inputfilename);
    ifstream infile( inputfilename.c_str() );

    if (infile.is_open())
    {
      string line;
      double x, y, TLocal, eLocal, ux, uy, pixx, piyy, pixy, pizz, pietaeta;
      while ( getline (infile, line) )
      {
        //if ( count++ < n_header_line ) continue;

        istringstream iss(line);
        iss >> x >> y >> TLocal >> ux >> uy >> pixx >> piyy >> pixy >> pizz;

        TLocal  /= hbarc_GeVfm;                           // 1/fm
        eLocal   = 3.0*cpLoc*TLocal*TLocal*TLocal*TLocal; // 1/fm^4
        pixx    /= hbarc_GeVfm;                           // 1/fm^4
        piyy    /= hbarc_GeVfm;                           // 1/fm^4
        pixy    /= hbarc_GeVfm;                           // 1/fm^4
        pietaeta = pizz/(tau0*tau0*hbarc_GeVfm);          // 1/fm^6

        // initialize particle object
        Particle p;
        p.r(0)           = x;
        p.r(1)           = y;
        p.input.e        = eLocal;
        p.hydro.u(0)     = ux;
        p.hydro.u(1)     = uy;
        p.hydro.shv(0,0) = 0.0;
        p.hydro.shv(0,1) = 0.0;
        p.hydro.shv(0,2) = 0.0;
        p.hydro.shv(1,0) = 0.0;
        p.hydro.shv(1,1) = pixx;
        p.hydro.shv(1,2) = pixy;
        p.hydro.shv(2,0) = 0.0;
        p.hydro.shv(2,1) = pixy;
        p.hydro.shv(2,2) = piyy;
        p.hydro.shv33    = pietaeta;

        systemPtr->particles.push_back( p );
      }

      infile.close();
    }
    else
    {
      std::cerr << "File " << inputfilename << " could not be opened!\n";
      abort();
    }
  }
  //----------------------------------------------------------------------------
  else
  {
      std::cerr << "Initial condition type " << initial_condition_type
                << " is not supported.\n";
      abort();
  }


  return;
}

//------------------------------------------------------------------------------
void InputOutput::print_system_state()
{
  //---------------------------------
  if (settingsPtr->printing_to_txt)
    print_system_state_to_txt();

  //---------------------------------
  if (settingsPtr->printing_to_HDF)
    print_system_state_to_HDF();

  //---------------------------------
  print_freeze_out();

  //---------------------------------
  // increment timestep index
  n_timesteps_output++;

  return;
}

//------------------------------------------------------------------------------
void InputOutput::print_system_state_to_txt()
{
  string outputfilename = output_directory + "/system_state_"
                          + std::to_string(n_timesteps_output) + ".dat";
  ofstream out( outputfilename.c_str() );

  out << systemPtr->t << "\n";
  int iParticle = 0;
  if ( settingsPtr->using_Gubser )
    for ( auto & p : systemPtr->particles )
    {
      out << p.r << " "
          << p.T() << " "
          << p.e() << " "
          << p.hydro.u(0) << " "
          << p.hydro.u(1) << " "
          << p.hydro.shv(1,1) << " "
          << p.hydro.shv(2,2) << " "
          << p.hydro.shv(1,2) << " "
          << pow(systemPtr->t,2.0)*p.hydro.shv33 << " "
          << p.rhoB() << " "
          << p.rhoS() << " "
          << p.rhoQ() << "\n";
      }
  else
  {
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
          /*<< p.hydro.bigPI << " "     //32
          << p.hydro.tauRelax << " "
          << p.hydro.Bulk << " "
          << p.hydro.dBulk_dt << " "
          << p.hydro.zeta << " "
          << p.hydro.dsigma_dt << " "
          << p.hydro.div_u << " "       //38
          << p.hydro.du_dt << "       "
          << p.hydro.gradV << "       "
          << p.hydro.gradU << "       "
          << p.hydro.gradBulk << "       "
          << p.hydro.gradshear << "       "
          << p.hydro.divshear << "   "
          << p.contribution_to_total_E << "   "
          << p.contribution_to_total_Ez << "   "*/
          << p.get_current_eos_name() << "\n";

  }

  out << std::flush;
  
  out.close();

  return;
}

//------------------------------------------------------------------------------
void InputOutput::print_system_state_to_HDF()
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
  vector<double> parameters      = { systemPtr->t, systemPtr->e_2_X.back(),
                                                   systemPtr->e_2_P.back() };

  hdf5_file.output_dataset( dataset_names, dataset_units, data, width,
                            n_timesteps_output, eos_tags,
                            parameters, parameter_names );

  return;
}





//==============================================================================
void InputOutput::print_freeze_out()
{
  string outputfilename = output_directory + "/freeze_out_"
                          + std::to_string(n_timesteps_output) + ".dat";
  ofstream FO( outputfilename.c_str(), ios::out | ios::app );

  auto & fo = wsPtr->fo;

  if ( fo.divTtemp.size() != fo.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( fo.divT.size() != fo.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( fo.gsub.size() != fo.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( fo.uout.size() != fo.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( fo.swsub.size() != fo.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( fo.bulksub.size() != fo.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( fo.shearsub.size() != fo.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( fo.shear33sub.size() != fo.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( fo.tlist.size() != fo.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( fo.rsub.size() != fo.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( fo.sFO.size() != fo.cf ) cout << "WARNING: WRONG SIZE" << endl;
  if ( fo.Tfluc.size() != fo.cf ) cout << "WARNING: WRONG SIZE" << endl;

  for (int i = 0; i < fo.cf; i++)
    FO << fo.divTtemp[i] << " " << fo.divT[i] << " "
        << fo.gsub[i] << " " << fo.uout[i] << " "
        << fo.swsub[i] << " " << fo.bulksub[i] << " " 
        << fo.shearsub[i](0,0) << " "
        << fo.shearsub[i](1,1) << " " 
        << fo.shearsub[i](2,2) << " "
        << fo.shear33sub[i] << " " 
        << fo.shearsub[i](1,2) << " " 
        << fo.tlist[i] << " " << fo.rsub[i] << " "
        << fo.sFO[i] << " " << fo.Tfluc[i] << endl;

  FO.close();

  return;
}
