#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>

//using namespace std;
using std::cout;
using std::endl;
using std::string;

#include "../include/constants.h"
#include "../include/eos.h"
#include "../include/formatted_output.h"
#include "../include/linklist.h"
#include "../include/particle.h"
#include "../include/stopwatch.h"
#include "../include/system_state.h"
#include "../include/vector.h"

using namespace constants;

////////////////////////////////////////////////////////////////////////////////
void SystemState::set_SettingsPtr(Settings * settingsPtr_in)
{
  settingsPtr = settingsPtr_in;
}

////////////////////////////////////////////////////////////////////////////////
void SystemState::initialize()  // formerly called "manualenter"
{
  formatted_output::report("Initializing system");

  n_particles = particles.size();
  t           = settingsPtr->t0;
  h           = settingsPtr->h;

  formatted_output::update("set freeze out parameters");

  formatted_output::detail("freeze out temperature = "
                           + to_string(settingsPtr->Freeze_Out_Temperature
                                        *hbarc_MeVfm) + " MeV");
  formatted_output::detail("freeze out energy density = "
                           + to_string(efcheck*hbarc_GeVfm) + " GeV/fm^3");
  formatted_output::detail("freeze out entropy density = "
                           + to_string(sfcheck) + " 1/fm^3");

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
  formatted_output::report("Initializing linklist");

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
    S += p.specific.s*p.norm_spec.s;

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
    Btotal += p.specific.rhoB*p.norm_spec.rhoB;
    Stotal += p.specific.rhoS*p.norm_spec.rhoS;
    Qtotal += p.specific.rhoQ*p.norm_spec.rhoQ;
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
           * p.norm_spec.s * t / p.hydro.sigma;
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
           * p.norm_spec.s / p.hydro.sigma;
    dEz += p.contribution_to_total_dEz;
  }

}



