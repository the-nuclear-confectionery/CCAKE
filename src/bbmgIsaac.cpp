#ifndef BBMG_CPP
#define BBMG_CPP

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include "../include/constants.h"
#include "../include/bbmgheader.h"
//#include "kernel.h"
//#include "particle.h"
//#include "settings.h"
//#include "system_state.h"

using namespace std;

BBMG::BBMG( Settings * settingsPtr_in, SystemState * systemPtr_in )
  : settingsPtr{ settingsPtr_in },
    systemPtr{ systemPtr_in }
{
  srand( time( NULL ) );
  TD    = 150;// Decoherence temperature in MeV
  q     = 1;// Fluctuation parameter, subject to change
  Cg    = 3; // Cassimir const gluons
  Cq    = 4./3; // Cassimir const quarks
  z     = 1; // path length dependence
  a     = 0; //Initial jet energy dependence
  c     = (2+z-a)/3; //medium temperature dependence
   // Need to use pointer to get T dependence or just use Matt's eqn with the access to ff[i].T
  vjet  = 1;
  area  = PI*pow(2.*systemPtr->h,2);
  rr.resize(systemPtr->n());

  for (int i = 0; i < 15; i++)
  {
    Rq[i]  = 0;
    Rg[i]  = 0;
    phi[i] = i*PI/7;
  }

  Pfg = 10;
  Pfq = 10;

  //call to complete initialization
  //initial();
}


double BBMG::get_kappa(double T) //This is taken from one plot we found of Barbara's, no reason to believe this should be set in stone. Will test multiple Kappa's
{ 
  return 2.5*exp(-9.22*T) + 0.04;
}

 
void BBMG::initial()
{ 
  rho0tot = 0;
  for ( int i = 0; i < systemPtr->particles.size(); ++i )
  {
    auto & p = systemPtr->particles[i];
    double rsub = p.p() / p.T();
    rho0tot    += rsub;
    //cout << "p.T is " << p.T()*constants::hbarc_MeVfm << " and TD is " << TD << endl;
    if ( p.T() * constants::hbarc_MeVfm > TD )
    {
      field sub; //field of all particles and their info, eventually turned into vector ff for calculations
      sub.r[0] = p.r(0);
      //cout << "Positions (x) of each particle in the grid is " << p.r(0) << "\n";
      sub.r[1] = p.r(1);
      sub.rho0 = rsub; //Density left in terms of femtometers
      sub.sph  = i;
      sub.T    = p.T() * constants::hbarc_MeVfm;
      cout << endl << "This is the value of sub.T for each particle: " << sub.T << endl;
      
      double kappa = get_kappa(p.T() * constants::hbarc_MeVfm);
      sub.line = 0.5 * kappa * pow(settingsPtr->t0, z) * pow(sub.rho0, c) * systemPtr->dt; // only if initial flow=0

      for (int j=0; j<14; j++) //initializes jets at each point in grid space, over 14 directions
      {
        sub.phi = phi[j];
        sub.pid = j;
        sub.on  = 1;
      }
      ff.push_back(sub);
    }
  }
}



double BBMG::flow(field &f) { return f.gam*(1-f.vmag*cos(f.phi-f.vang)); }


double BBMG::gft(double p) { return 2*p; } 


double BBMG::qft(double p) { return 2*p; }


double BBMG::efluc()
{
  int random_variable = std::rand()/RAND_MAX;
  double zeta         = random_variable*(q+2.);
  return (1.+q) / pow(q+2, 1+q) * pow(q+2.-zeta, q);
}


void BBMG::propagate()
{
  double tau  = systemPtr->t + settingsPtr->t0;
  int stillon = 0;
  int tot     = ff.size();
  double P0g = 0, P0q = 0;
  for (int i = 0; i < tot; i++)
  {
    // propagate x,y position of jet
    ff[i].r[0] += vjet * systemPtr->dt * cos(ff[i].phi);
    ff[i].r[1] += vjet * systemPtr->dt * sin(ff[i].phi);

    double kappa = get_kappa(ff[i].T);
      
    inter( ff[i] ); //interpolation of the field
    cout << "ff[i].T is: " << ff[i].T << "\n"; //for some reason, still higher than expected. still checking how interpolation is working
    if ( ( ff[i].on == 1 ) && ( ff[i].T > TD ) )
    {
      ff[i].line += pow(tau, z) * pow(ff[i].rho, c) * systemPtr->dt; // * flow(ff[i])
      stillon++;
    }
    else
    {
      ff[i].on    = 0;
      ff[i].line += 0.5 * kappa * pow(tau,z) * pow(ff[i].rho, c) * systemPtr->dt; /* flow(ff[i])*/
      //ff[i].line *= efluc();

      P0g  = Pfg + Cg * ff[i].line; //* pow(Pfg, 1-a)
      P0q  = Pfq + Cq * ff[i].line; //* pow(Pfq, 1-a) 
      cout << "This is the value of P0g: " << P0g << "This is the value of P0q: " << P0q << endl;

      int jj      = ff[i].pid;
      //Rq[jj]     += pow(P0g/Pfg, 1+a) * ff[i].rho0 * gft(P0g) / gft(Pfg);
      //Rq[jj]     += pow(P0q/Pfq, 1+a) * ff[i].rho0 * qft(P0g) / qft(Pfg); 
    }
  }

  if ( stillon == 0 )
  {
    for (int j=0; j<14; j++)
    {
      Rq[j] /= rho0tot;
      Rg[j] /= rho0tot;
    }
  }
}


void BBMG::inter( field &f ) //How are f.T and others working here? It seems to be rewriting sub.T from above?
{
  /*if (initial_condition_type == "ICCING")
  {
  }
  else
  {
  }*/
  double den = 0, den2 = 0;
  for ( auto & p : systemPtr->particles )
  {
    double dx    = p.r(0)-f.r[0];
    //cout << "X position from f.r is " << f.r[0] << "\n";
    double dy    = p.r(1)-f.r[1];

    double rdiff = sqrt(dx*dx+dy*dy)/systemPtr->h;

    if (rdiff<2)
    {
      ++den;
      den2     += p.norm_spec.s;
      double kk = kernel::kernel(rdiff); //In order to call this, shouldn't I be including kernel.cpp instead of kernel.h?"
      double gridx = settingsPtr->stepx;
      double gridy = settingsPtr->stepy;
      //f.T      += p.T()*constants::hbarc_MeVfm*0.06*0.06*kk; //Here I changed the hardcoded grid size to a read in of the default grid from settings.
      f.T      += p.T()*constants::hbarc_MeVfm*gridx*gridy*kk; //Interpolation seemingly done in fm instead of MeV? After correcting with the constants list, almost correct
      //cout << "Interpolated temp value is " << f.T << "\n";
      f.rho    += (p.p()/p.T())*kk;
      f.v[0]   += p.hydro.v(0)*kk;
      f.v[1]   += p.hydro.v(1)*kk;

      //cout << dx << " " << dy << " " << p.T()*constants::hbarc_MeVfm << " " << p.hydro.v << endl;
    }
  }

  double fac = den / area;
  //f.T       *= constants::hbarc_MeVfm;
  f.rho     /= fac;
  f.v[0]    /= fac;
  f.v[1]    /= fac;
  f.vmag     = sqrt( f.v[0]*f.v[0] + f.v[1]*f.v[1] );
  f.vang     = atan2( f.v[1], f.v[0] );
  f.gam      = 1.0 / sqrt( f.vmag*f.vmag + 1.0 );
}


#endif
