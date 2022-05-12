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
#include "eos.h"
#include "linklist.h"
#include "particle.h"
#include "Stopwatch.h"
#include "system_state.h"
#include "vector.h"

using namespace constants;


////////////////////////////////////////////////////////////////////////////////
void SystemState::set_SettingsPtr(Settings * settingsPtr_in)
{
  settingsPtr = settingsPtr_in;
}


////////////////////////////////////////////////////////////////////////////////
void SystemState::initialize()  // formerly called "manualenter"
{
  n_particles = particles.size();
  t           = settingsPtr->t0;
  h           = settingsPtr->h;

  std::cout << "FO temp. = " << settingsPtr->Freeze_Out_Temperature << " 1/fm\n";
  std::cout << "efcheck = " << efcheck*hbarc_GeVfm << " GeV/fm^3\n";
  std::cout << "sfcheck = " << sfcheck << " 1/fm^3\n";

  // initialize information for particles
  for (auto & p : particles)
  {
    p.set_SettingsPtr( settingsPtr );
    p.efcheck = efcheck;
  }

  return;
}



void SystemState::initialize_linklist()
{
  cout << "Initial conditions type: " << settingsPtr->IC_type << endl;

  // initialize linklist
  linklist.initialize( &particles, settingsPtr->h );

  // this should be moved
  // use a lambda function for convenience
  for (auto & p : particles)
    p.hydro.g2 = [](double x){ return x*x; }( p.gamcalc() );

  return;
}









///////////////////////////////////////
void SystemState::conservation_entropy()
{
  S = 0.0;
  for ( auto & p : particles )
    S += p.specific.s*p.sigmaweight;

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
    Btotal += p.smoothed.rhoB*p.rhoB_weight;
    Stotal += p.smoothed.rhoS*p.rhoS_weight;
    Qtotal += p.smoothed.rhoQ*p.rhoQ_weight;
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
         = ( p.hydro.C*p.hydro.g2 - p.p() - p.hydro.bigPI + p.hydro.shv(0,0) )
           * p.sigmaweight * t / p.hydro.sigma;
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
         = ( p.p() + p.hydro.bigPI + p.hydro.shv33*t2 )
           * p.sigmaweight / p.hydro.sigma;
    dEz += p.contribution_to_total_dEz;
  }

}



//==============================================================================
// this routine is used to initialize quantities prior to RK evolution
void SystemState::set_current_timestep_quantities()
{
  etasigma0.resize(n_particles);
  Bulk0.resize(n_particles);
  particles_E0.resize(n_particles);

  u0.resize(n_particles);
  r0.resize(n_particles);

  shv0.resize(n_particles);

  for (int i = 0; i < n_particles; ++i)
  {
    auto & p = particles[i];

    u0[i]        = p.hydro.u;
    r0[i]        = p.r;
    etasigma0[i] = p.specific.s;
    Bulk0[i]     = p.hydro.Bulk;
    mini( shv0[i], p.hydro.shv );

    particles_E0[i] = p.contribution_to_total_Ez;
  }
}

