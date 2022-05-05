#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <sstream>

//using namespace std;
using std::cout;
using std::endl;
using std::string;
using std::stod;

#include "constants.h"
#include "vector.h"
#include "particle.h"
#include "eos.h"
#include "Stopwatch.h"
#include "system_state.h"
#include "linklist.h"

using namespace constants;


////////////////////////////////////////////////////////////////////////////////
void SystemState::set_SettingsPtr(Settings * settingsPtr_in)
{
  settingsPtr = settingsPtr_in;
}


////////////////////////////////////////////////////////////////////////////////
void SystemState::initialize()  // formerly called "manualenter"
{
  int start, end;
  int df;

  cout << __PRETTY_FUNCTION__ << ": " << endl;
  _n = particles.size();

  cout << "_n = " << _n << endl;

  t = settingsPtr->t0;

  cout << "t = " << t << endl;

  _h = settingsPtr->_h;

  // set viscosities (this will have to change when bringing in transport
  // coefficients class)
  // (probably add if statement to check if they're constants or variable)
  std::cout << "settingsPtr->etaOption = " << settingsPtr->etaOption << std::endl;
  std::cout << "settingsPtr->zetaOption = " << settingsPtr->zetaOption << std::endl;
  svf = stod(settingsPtr->etaOption);
  bvf = stod(settingsPtr->zetaOption);

  std::cout << "FO temp. = " << settingsPtr->Freeze_Out_Temperature << " 1/fm\n";
  std::cout << "efcheck = " << efcheck*hbarc_GeVfm << " GeV/fm^3\n";
  std::cout << "sfcheck = " << sfcheck << " 1/fm^3\n";

  freezeoutT = settingsPtr->Freeze_Out_Temperature;


  // initialize information for particles
  for (auto & p : particles)
  {
    p.set_SettingsPtr( settingsPtr );
    p.freezeoutT = freezeoutT;
    p.efcheck = efcheck;
  }

  linklist.efcheck = efcheck;
  linklist.sfcheck = sfcheck;
  linklist.fcount  = 0;
  linklist.average = 0;

  return;
}



void SystemState::initialize_linklist()
{
  cout << "Initial conditions type: " << settingsPtr->IC_type << endl;

  if ( settingsPtr->IC_type == "ICCING" )
  {
    settingsPtr->gtyp   = 6;

    int count           = 1;
    vector<string>        filelist( count );

    int j               = 0;
    filelist[j]         = "./ic0.dat"; // only doing single event
    linklist.filenames  = filelist;
    linklist.fcount     = count;
    linklist.fnum       = linklist.start;
    
    int currently_frozen_out = number_part;
    linklist.initialize( settingsPtr->t0, particles.size(),
                         settingsPtr->_h, &particles, dt, currently_frozen_out );

    linklist.gtyp=settingsPtr->gtyp;

  }
  else if (    settingsPtr->IC_type == "Gubser"
            || settingsPtr->IC_type == "Gubser_with_shear" )
  {
    settingsPtr->gtyp   = 7;

    int count           = 1;
    vector<string>        filelist( count );

    int j               = 0;
    filelist[j]         = "./ic0.dat"; // only doing single event
    linklist.filenames  = filelist;
    linklist.fcount     = count;
    linklist.fnum       = linklist.start;
    
    int currently_frozen_out = number_part;
    linklist.initialize( settingsPtr->t0, particles.size(),
                         settingsPtr->_h, &particles, dt, currently_frozen_out );

    linklist.gtyp=settingsPtr->gtyp;
  }
  else
  {
    std::cerr << "Initial conditions type = " << settingsPtr->IC_type
              << " not supported!" << std::endl;
    exit(1);
  }


  for (auto & p : particles)
  {
    double gg = p.gamcalc();
    p.g2      = gg*gg;
  }

  return;
}









///////////////////////////////////////
void SystemState::conservation_entropy()
{
  S = 0.0;
  for ( auto & p : particles )
    S += p.eta_sigma*p.sigmaweight;

  if (linklist.first==1) S0 = S;
}

///////////////////////////////////////
void SystemState::conservation_BSQ()
{
  // reset
  Btotal = 0.0;
  Stotal = 0.0;
  Qtotal = 0.0;

  // sum
  for ( auto & p : particles )
  {
    Btotal += p.rhoB_sub*p.rhoB_weight;
    Stotal += p.rhoS_sub*p.rhoS_weight;
    Qtotal += p.rhoQ_sub*p.rhoQ_weight;
  }

  // save initial totals
  if (linklist.first==1)
  {
    Btotal0 = Btotal;
    Stotal0 = Stotal;
    Qtotal0 = Qtotal;
  }
  return;
}




///////////////////////////////////////
void SystemState::conservation_energy()
{
  ///////////////////////////////////////////////
  // don't bother checking energy conservation on
  // intermediate RK steps
  if ( rk2 == 1 )
  {
    // calculate total energy (T^{00})
    E = 0.0;
    for ( auto & p : particles )
    {
      p.contribution_to_total_E
         = ( p.C*p.g2 - p.p() - p.bigPI + p.shv(0,0) )
           * p.sigmaweight * t / p.sigma;
      E += p.contribution_to_total_E;
    }

    // store initial total energy
    // for checking subsequent energy loss
    if (linklist.first==1)
    {
      linklist.first = 0;
      E0             = E;
    }

    // Ez is initially set to zero,
    // updated subsequently during RK integration
    Etot  = E + Ez;
    Eloss = (E0-Etot)/E0*100;
    rk2   = 0;
  }

  ///////////////////////////////////////////////
  // this enters the RK routine and should be
  // done for intermediate steps as well;
  // this gives the longitudinal energy flux (~T^{\eta\eta})
  dEz = 0.0;
  double t2 = t*t;
  for ( auto & p : particles )
  {
    p.contribution_to_total_dEz
         = ( p.p() + p.bigPI + p.shv33*t2 )
           * p.sigmaweight / p.sigma;
    dEz += p.contribution_to_total_dEz;
  }

}



///////////////////////////////////////////////////////////////////////////////
void SystemState::set_current_timestep_quantities()
{
  N = _n;

  etasigma0.resize(N);
  Bulk0.resize(N);
  particles_E0.resize(N);

  u0.resize(N);
  r0.resize(N);

  shv0.resize(N);

  for (int i=0; i<N; ++i)
  {
    auto & p = particles[i];

    u0[i]        = p.u;
    r0[i]        = p.r;
    etasigma0[i] = p.eta_sigma;
    Bulk0[i]     = p.Bulk;
    mini( shv0[i], p.shv );

    particles_E0[i] = p.contribution_to_total_Ez;
  }
}


///////////////////////////////////////////////////////////////////////////////
void SystemState::get_derivative_step(double dx)
{
  N = _n;

  for (int i=0; i<N; ++i)
  {
    auto & p = particles[i];

    p.r            = r0[i]        + dx*p.v;

    if ( p.Freeze < 5 )
    {
      p.u            = u0[i]        + dx*p.du_dt;
      p.eta_sigma    = etasigma0[i] + dx*p.detasigma_dt;
      p.Bulk         = Bulk0[i]     + dx*p.dBulk_dt;
      tmini( p.shv,    shv0[i]      + dx*p.dshv_dt );
    }
  }
}


