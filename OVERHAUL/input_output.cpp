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

#include "mathdef.h"
#include "vector.h"
#include "particle.h"
#include "input_output.h"
#include "constants.h"

using namespace constants;

using std::endl;
using std::flush;
using std::string;

// Constructors and destructors.
InputOutput::InputOutput(){}
InputOutput::~InputOutput(){}

void InputOutput::set_EquationOfStatePtr( EquationOfState * eosPtr_in )
{
  eosPtr = eosPtr_in;
}

void InputOutput::set_EquationsOfMotionPtr( EquationsOfMotion * eomPtr_in )
{
  eomPtr = eomPtr_in;
}

void InputOutput::set_SettingsPtr( Settings * settingsPtr_in )
{
  settingsPtr = settingsPtr_in;
}

void InputOutput::set_SystemStatePtr( SystemState * systemPtr_in )
{
  systemPtr = systemPtr_in;
}


void InputOutput::set_results_directory( string path_to_results_directory )
{
  output_directory = path_to_results_directory;
}



void InputOutput::load_settings_file( string path_to_settings_file )
{
    string Param_file = path_to_settings_file;
    ifstream infile( Param_file.c_str() );
    if (infile.is_open())
    {
        string line;
        string ignore = "";
        string param = "";
        vector<string> all_parameters;
        while ( getline (infile, line) )
        {
            istringstream iss(line);
            iss >> ignore >> param;
            all_parameters.push_back(param);
        }

cout << "all_parameters.size() = " << all_parameters.size() << endl;
for ( auto & entry : all_parameters )
  cout << entry << endl;

        settingsPtr->IC_type                = all_parameters[0];
        settingsPtr->_h                     = stod(all_parameters[1]);
        settingsPtr->dt                     = stod(all_parameters[2]);
        settingsPtr->t0                     = stod(all_parameters[3]);
        settingsPtr->EoS_type               = all_parameters[4];
        settingsPtr->EoS_option             = all_parameters[5];
        settingsPtr->eta                    = all_parameters[6];
        settingsPtr->zeta                   = all_parameters[7];
        settingsPtr->Freeze_Out_Temperature = stod(all_parameters[8])/hbarc_MeVfm;
        settingsPtr->Freeze_Out_Type        = all_parameters[9];

        infile.close();
    }

    return;
}

void InputOutput::set_EoS_type()
{
  string EoS_type = settingsPtr->EoS_type;
  string EoS_option = settingsPtr->EoS_option;
  string EoS_files_location = "EoS/" + EoS_type + "/" + EoS_option;
  string densities = EoS_files_location + "/densities.dat";
  string derivatives = EoS_files_location + "/derivatives.dat";

  if (EoS_option == "Default")
  {
    cout << "Running default EoS option for " << EoS_type << endl;
  }
  else
  {
    cout <<"EoS option not recognized for " << EoS_type << ", now exiting." << endl;
    exit(1);
  }


  eosPtr->quantity_file = densities;
  eosPtr->deriv_file = derivatives;
  return;
}

void InputOutput::read_in_initial_conditions()
{
  string initial_condition_type = settingsPtr->IC_type;
  int total_header_lines;
  string IC_file = "initial_conditions/";

  if (initial_condition_type == "ICCING")
  {
      cout << "Reading in ICCING initial conditions!" << endl;
      IC_file = IC_file+"/Iccing_conditions.dat"; // need to change ic0.dat
      total_header_lines = 1;    
  }
  else
  {
      cout << "Selected initial condition type not supported." << endl;
      exit(1);
  }

  ifstream infile(IC_file.c_str());
  cout << "Initial conditions file: " << IC_file << endl;
  if (infile.is_open())
  {
    string line;
    int count_header_lines = 0;
    int count_file_lines   = 0;
    double x, y, e, rhoB, rhoS, rhoQ;
    double ignore, stepX, stepY;

    while (getline (infile, line))
    {
      istringstream iss(line);
      if(count_header_lines < total_header_lines)
      {
        settingsPtr->headers.push_back(line);
        iss >> ignore >> stepX >> stepY;
        settingsPtr->stepx = stepX;
        settingsPtr->stepy = stepY;
        count_header_lines++;
      }
      else
      {
        iss >> x >> y >> e >> rhoB >> rhoS >> rhoQ;
        e /= hbarc_GeVfm;
        vector<double> fields({x,y,e,rhoB,rhoS,rhoQ});
        systemPtr->particles.push_back( Particle(fields) );
      }
    }
  }
  else
  {

    cout << "Can't open " << IC_file << endl;
    exit(1);
  }
  infile.close();
}








void InputOutput::print_system_state()
{
  string outputfilename = output_directory + "/system_state_"
                          + std::to_string(n_timesteps_output) + ".dat";
  ofstream out( outputfilename.c_str() );

  out << systemPtr->t << endl;
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
				<< p.e()*hbarc_MeVfm << " "
				<< p.rhoB() << " "
				<< p.rhoS() << " "
				<< p.rhoQ() << " "
				<< p.s() << " "
				<< p.eta/(p.gamma*systemPtr->t) << " "
				<< p.eta_sigma << " "
				<< p.sigma << " " 
				<< p.sigmaweight << " "
				<< p.stauRelax << " " 
				<< p.bigtheta << " "
				<< sqrt( p.shv.x[0][0]*p.shv.x[0][0]
                -2.0*p.shv.x[0][1]*p.shv.x[0][1]
                -2.0*p.shv.x[0][2]*p.shv.x[0][2]
                + p.shv.x[1][1]*p.shv.x[1][1]
                + p.shv.x[2][2]*p.shv.x[2][2]
                +2.0*p.shv.x[1][2]*p.shv.x[1][2]
                +pow(systemPtr->t,4.0)*p.shv33*p.shv33 ) << " "
				<< p.stauRelax/systemPtr->t * p.bigtheta << " "
				<< p.u.x[0]/p.gamma << " "
				<< p.u.x[1]/p.gamma << " "
				<< p.gamma << endl;
  
  out.close();

  // increment timestep index
  n_timesteps_output++;

  return;
}