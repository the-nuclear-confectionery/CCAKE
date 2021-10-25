#ifndef _ENTERIC_CPP_
#define _ENTERIC_CPP_

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <vector>
#include "mathdef.h"
#include "vector.h"
#include "tables.h"
#include "particle.h"
#include "enteric.h"
#include<dirent.h>
#include <stdlib.h>
#include <algorithm>
#include <cmath>
#include "input_output.h"

// Constructors and destructors.
  input_output::input_output(){}
 ~input_output::input_output(){}


void input_output::read_in_initial_conditions()
{
 
    
  
}

void input_output::load_settings_file( string path_to_settings_file )
{
    string Param_file = path_to_settings_file+"Input_Parameters.inp";
    ifstream infile( Param_file.c_str() );
    if (infile.is_open())
    {
        string line;
        string ignore, param;
        vector<string> all_parameters;
        while ( getline (infile, line) )
        {
            istringstream iss(line);
            iss >> ignore >> param;
            all_parameters.push_back(param)
        }

        input_parameters.IC_type                = all_parameters[0];
        input_parameters.h                      = stod(all_parameters[1]);
        input_parameters.dt                     = stod(all_parameters[2]);
        input_parameters.t0                     = stod(all_parameters[3]);
        input_parameters.EoS_type               = param[4];
        input_parameters.EoS_option             = param[5];
        input_parameters.eta                    = param[6]
        input_parameters.zeta                   = param[7]
        input_parameters.Freeze_Out_Temperature = stod(param[8])
        input_parameters.Freeze_Out_Type        = param[9]

        infile.close();
    }
    return
}

void input_output::set_EoS_type()
{
  EoS_type = input_parameters.EoS_type;
  string EoS_files_location = 'EoS/' + EoS_type;
  string densities = EoS_files_location + '/densities.dat';
  string derivatives = EoS_files_location + '/derivatives.dat';
  string EoS_option = input_paramters.EoS_option;

  switch(EoS_spec)
  {
    Case default   :
      cout << "Running default EoS option for " << EoS_type << endl;
  }

  eos.innit(desnities,derivatives);
  return
}

void input_output::read_in_initial_conditions()
{
  string initial_condition_type = input_parameters.IC_type;
  int n_header_lines;
  string IC_file = 'All_Initial_Conds/'
  switch(initial_condition_type)
  {
    case 'ICCING' :
      cout << "Reading in ICCING initial conditions!" << endl;
      IC_file = IC_file+'Iccing_conditions.dat' // need to change ic0.dat
    case default :
      cout << "Selected initial condition type not supported."
  }

ifstream infile(IC_file.c_str());
cout << "Initial conditions file: " << IC_file << endl;
if (infile.is_open())
{
  string line;
  vector<vector<double>> all_initials;
  int count=0;
  if count<n_header_lines
while ( getline (infile, line) )
        {
            istringstream iss(line);
            iss >> xsub >> ysub >> esub >> rBsub >> rSsub >> rQsub;

        }
else
{
 if (y[2]>0.00301)
//if (y[2]>max(0.01,0.5*efcheck*hbarC)) //N.B. leave wiggle room around FO temp
{
// check if this file has finite energy density and zero charge densities
// (seems to cause problems in EOS)
const double eps_local = 1e-6;
/*if (   y[2] >= eps_local && y[3] < eps_local
	&& y[4] < eps_local  && y[5] < eps_local)
	{
		 y[3] = eps_local;
		 y[4] = eps_local;
		 y[5] = eps_local;
	}*/
iss >> var1;
    xsub.push_back(y[0]);
    iss>>var2;
    ysub.push_back(y[1]);
    iss>>var3;
    esub.push_back(y[2]);
    iss>>var4;
    rBsub.push_back(y[3]);
    iss>>var5;
    rSsub.push_back(y[4]);
    iss>>var6;
    rQsub.push_back(y[5]);
}
count ++;

if (!infile.is_open())
    {
        cout << "Can't open " << IC_file << endl;
        exit(1);
    }
  }
}
}
{
//cout << "CHECK(" << __LINE__ << "): " << j << "   " << x[j] << "   " << y[j] << endl;
        }
}
infile.close();

// Here goes the formation of SPH from the initial densities of ICCING 
numberof_sph=xsub.size();
p= new Particle<2>[numberof_sph];

cout << "After e-cutoff and freeze-out: size = " << numberof_sph << endl;
int kk=numberof_sph;
initial_frozen_sph=0;	//number of frozen out particles initially
for(int j=0; j<numberof_sph; j++)
	{
const auto & p = _p[j]; 
p.r.x[0]=xsub[j]; // set the initial x position of an sph particle 
p.r.x[1]=ysub[j]; // set the initial y position of an sph particle 
// p.e_sub=EOS.e_out(factor*esub[j]);
p.e_sub=esub[j]/hbarC; // set the initial energy density of an sph particle       
// not s_an!!  convert to 1/fm^4 and do not rescale by factor!
p.u.x[0]=0; // set the initial flow velocity
p.u.x[1]=0; // 
p.entropy_sigma = 1;// it was eta_sigma
p.sigmaweight=stepx*stepy; // nu_\alpha
p.rho_weight = stepx*stepy; // Chris?
p.Bulk = 0;
p.B=rBsub[j]*stepx*stepy;	// nu alpha * 		
p.S=rSsub[j]*stepx*stepy;			
p.Q=rQsub[j]*stepx*stepy;			
p.rhoB_an=rBsub[j];					
p.rhoS_an=rSsub[j];					
p.rhoQ_an=rQsub[j];					
p.transverse_area = stepx*stepy;

if (j==0)
	cout << "readICs_iccing(" << __LINE__ << "): "
	<< "SPH particles: "
	<< j << "   "
	<< p.r.x[0] << "   " << p.r.x[1] << "   "
	<< p.e_sub << "   " << p.rhoB_an << "   "
	<< p.rhoS_an << "   " << p.rhoQ_an << "   "
	<< p.sigmaweight << endl;
// make educated initial guess here for this particle's (T, mu_i) coordinates
// (improve this in the future)
	p.SPH_cell.T   = 500.0/hbarc_MeVfm;	// rootfinder seems to work better going downhill than "uphill"
	p.SPH_cell.muB = 0.0/hbarc_MeVfm;
	p.SPH_cell.muS = 0.0/hbarc_MeVfm;
	p.SPH_cell.muQ = 0.0/hbarc_MeVfm;

  if (p.e_sub>efcheck)	          // impose freeze-out check for e, not s
  {
            p.Freeze=0;
        }
        else
        {
            p.Freeze=4;
            --kk;
            ++initial_frozen_sph;
        }
    }
cout << "After freezeout (redundant): size = " << numberof_sph-initial_frozen_sph << endl;
}
}

#endif