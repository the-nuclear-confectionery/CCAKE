#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <vector>
#include "mathdef.h"
#include "vector.h"
#include "particle.h"
#include <stdlib.h>
#include <algorithm>
#include <cmath>
#include "input_output.h"
#include "constants.h"

using namespace constants;

using std::flush;

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

        settingsPtr->input_parameters.IC_type                = all_parameters[0];
        settingsPtr->input_parameters.h                      = stod(all_parameters[1]);
        settingsPtr->input_parameters.dt                     = stod(all_parameters[2]);
        settingsPtr->input_parameters.t0                     = stod(all_parameters[3]);
        settingsPtr->input_parameters.EoS_type               = all_parameters[4];
        settingsPtr->input_parameters.EoS_option             = all_parameters[5];
        settingsPtr->input_parameters.eta                    = all_parameters[6];
        settingsPtr->input_parameters.zeta                   = all_parameters[7];
        settingsPtr->input_parameters.Freeze_Out_Temperature = stod(all_parameters[8]);
        settingsPtr->input_parameters.Freeze_Out_Type        = all_parameters[9];

        infile.close();
    }

    return;
}

void InputOutput::set_EoS_type()
{
  string EoS_type = settingsPtr->input_parameters.EoS_type;
  string EoS_files_location = "EoS/" + EoS_type;
  string densities = EoS_files_location + "/densities.dat";
  string derivatives = EoS_files_location + "/derivatives.dat";
  string EoS_option = settingsPtr->input_parameters.EoS_option;

  if(EoS_option == "Default")
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
  string initial_condition_type = settingsPtr->input_parameters.IC_type;
  int total_header_lines;
  string IC_file = "initial_conditions/";

  if(initial_condition_type == "ICCING")
  {
      cout << "Reading in ICCING initial conditions!" << endl;
      IC_file = IC_file+"Iccing_conditions.dat"; // need to change ic0.dat
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
    int count_file_lines = 0;
    double x,y,e,rhoB,rhoS,rhoQ;
    double ignore, stepX, stepY;
    istringstream iss(line);
    while (getline (infile, line))
            {
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
          auto & test = systemPtr->particles;
for ( auto & p : systemPtr->particles )
  cout << "new check: " << flush
      << flush << p.e_sub << "   "
      << flush << p.rhoB_an << "   "
      << flush << p.rhoS_an << "   "
      << flush << p.rhoQ_an << "   " << endl;
          cout << "new check: " << flush << test.size() << "   "
                << flush << (test[0]).e_sub << "   "
                << flush << (test[0]).rhoB_an << "   "
                << flush << (test[0]).rhoS_an << "   "
                << flush << (test[0]).rhoQ_an << "   " << flush;
          cout << (systemPtr->particles)[0].r.x[0] << "   "
                << (systemPtr->particles)[0].r.x[1] << endl;
  }
  else
  {

    cout << "Can't open " << IC_file << endl;
    exit(1);
  }
  infile.close();
}
