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
  Freezeout_Temp    = 150;// Decoherence temperature in MeV
  Cg                = 3; // Cassimir const gluons
  Cq                = 4./3; // Cassimir const quarks
  // All next quantities are part of the BBMG parameters
  q                 = 1;// Fluctuation parameter, subject to change
  z                 = 1; // path length dependence
  a                 = 0; //Initial jet energy dependence
  c                 = (2+z-a)/3; //medium temperature dependence
  
  //===============================================
  vjet  = 1;
  // area is taken from parameter h read in from input parameters file
  area  = PI*pow(2.*systemPtr->h,2);
  rr.resize(systemPtr->n());

  for (int i = 0; i < 15; i++)
  { // Not sure what this for yet but the quantities aren't even calculated
    Rq[i]  = 0;
    Rg[i]  = 0;
    phi[i] = i*PI/7;
  }
  //Setting final energy as a start point for the integration
  Pfg = 10000;
  Pfq = 10000;

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
    //Density from pressure over temperature
    double rsub = p.p() / p.T();
    //cout << "Pressures from particle list: " << p.p() << endl << "Temperatures from particle list: " << p.T() << endl;
    rho0tot    += rsub;
    if ( p.T() * constants::hbarc_MeVfm > Freezeout_Temp )
    {
      field sph_particle; //field of all sph particles where we take necessary line integral info
      sph_particle.r[0] = p.r(0);
      //cout << "Positions (x) of each particle in the grid is " << p.r(0) << "\n";
      sph_particle.r[1] = p.r(1);
      sph_particle.rho0 = rsub; //Density left in terms of femtometers
      sph_particle.sph  = i;
      sph_particle.T    = p.T() * constants::hbarc_MeVfm;
      //cout << endl << "Initial temps of non frozen-out sph particles in MeV is " << sph_particle.T << endl; This is working just fine
      
      double kappa = get_kappa(sph_particle.T);
      sph_particle.line = 0.5 * kappa * pow(settingsPtr->t0, z) * pow(sph_particle.rho0, c) * systemPtr->dt; // only if initial flow=0

      for (int j=0; j<14; j++) //initializes jets at each point in grid space, over 14 directions
      {
        sph_particle.phi = phi[j];
        sph_particle.pid = j;
        sph_particle.on  = 1;
      }
      full_sph_field.push_back(sph_particle);
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

//Delaunay triangulation if needed
/*void BBMG::initialize_Delaunay()
{

}*/

void BBMG::propagate()
{
  double tau  = systemPtr->t + settingsPtr->t0;
  int stillon = 0;
  int tot     = full_sph_field.size();
  double P0g = 0, P0q = 0;
  for (int i = 0; i < tot; i++)
  {
    // propagate x,y position of jet on top of sph particles
    full_sph_field[i].r[0] += vjet * systemPtr->dt * cos(full_sph_field[i].phi);
    full_sph_field[i].r[1] += vjet * systemPtr->dt * sin(full_sph_field[i].phi);


    double kappa = get_kappa(full_sph_field[i].T);
    //cout << "This is the value of our kappa coupling: " << kappa << endl;
    cout << "This is checking if tau is working properly: " << tau << endl;
    cout << "This is checking if the density is coming out positive: " << full_sph_field[i].rho0 << endl;
    inter( full_sph_field[i] ); //interpolation of the field
    //cout << "Interpolated field temp is: " << full_sph_field[i].T << "\n";
    //abort();

    
    if ( ( full_sph_field[i].on == 1 ) && ( full_sph_field[i].T > Freezeout_Temp ) )
    {
      full_sph_field[i].line += pow(tau, z) * pow(full_sph_field[i].rho0, c) * systemPtr->dt; // * flow(ff[i])
      cout << "Checking the values of the line integration: " << full_sph_field[i].line << endl;
      stillon++;
    }
    else //This comes in when we drop below freezeout temp, as .on should never go to 0 on its own
    {
      full_sph_field[i].on    = 0;
      full_sph_field[i].line += 0.5 * kappa * pow(tau,z) * pow(full_sph_field[i].rho0, c) * systemPtr->dt; /* flow(ff[i])*/
      //ff[i].line *= efluc();
      cout << "Checking values of line integration AFTER the jets hit FO Temperature: " << full_sph_field[i].line << endl;
      abort();
      P0g  = Pfg + Cg * full_sph_field[i].line; //* pow(Pfg, 1-a)
      P0q  = Pfq + Cq * full_sph_field[i].line; //* pow(Pfq, 1-a) 
      cout << "P0g: " << P0g << "MeV, P0q: " << P0q << "MeV" << endl;
      if ( P0g > 10000 || P0q > 10000 ) 
      {
        cout << "This is the value of P0g: " << P0g << "This is the value of P0q: " << P0q << endl;
        abort();
      }
      int jj      = full_sph_field[i].pid;
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


void BBMG::inter( field &f ) 
{
  /*if (initial_condition_type == "ICCING")
  {
  }
  else
  {
  }*/
  double norm = 0;
  double den = 0, den2 = 0;
  f.T    = 0;
  f.rho  = 0;
  f.v[0] = 0;
  f.v[1] = 0;
  for ( auto & p : systemPtr->particles )
  {
    double dx    = p.r(0)-f.r[0];
    //cout << "X position from f.r is " << f.r[0] << "\n" << "X position from p.r is " << p.r(0) << "\n";
    //cout << "dx is: " << dx << "\n";
    double dy    = p.r(1)-f.r[1];
    //cout << "Y position from f.r is " << f.r[1] << "\n" << "Y position from p.r is " << p.r(1) << "\n";
    //cout << "dy is: " << dy << "\n";

    double rdiff = sqrt(dx*dx+dy*dy)/systemPtr->h;
    //cout << "Distance used in kernel interpolation is: " << rdiff << "\n";  Returning values WAY above 2
    

    if (rdiff<2)
    {
      //cout << "Distance used in kernel interpolation is: " << rdiff << "\n";
      den++;
      den2     += p.norm_spec.s;
      double kern = kernel::kernel(rdiff); //In order to call this, shouldn't I be including kernel.cpp instead of kernel.h?"
      norm         += kern;
      double gridx = settingsPtr->stepx;
      double gridy = settingsPtr->stepy;
      //f.T      += p.T()*constants::hbarc_MeVfm*0.06*0.06*kk; //Here I changed the hardcoded grid size to a read in of the default grid from settings.
      f.T      += p.T()*constants::hbarc_MeVfm*kern; // After correcting with the constants list, almost correct --------- WHY IS THIS += AND NOT JUST =
      //cout << "Interpolated temp value is " << f.T << "\n";
      /*if (f.T > 900)
      {abort();}*/
      f.rho    += (p.p()/p.T())*kern;
      f.v[0]   += p.hydro.v(0)*kern;
      //cout << "X velocity is: " << f.v[0] << endl; This is fine for now, not seemingly above 1
      f.v[1]   += p.hydro.v(1)*kern;
      //By following the style of sph workstation, I have included den2 as the normalization factor with the kernel function


      //cout << dx << " " << dy << " " << p.T()*constants::hbarc_MeVfm << " " << p.hydro.v << endl;
      
    }
  }
  //This fac quantity seems to be fools gold, every time I have used it I get either 0's or no values at all
 //may need to move these inside the for loop?? including these returns 0 for all quantities below
  //f.T       *= constants::hbarc_MeVfm;
  f.T       /= norm;
  f.rho     /= norm;
  f.v[0]    /= norm;
  f.v[1]    /= norm;
  f.vmag     = sqrt( f.v[0]*f.v[0] + f.v[1]*f.v[1] );
  f.vang     = atan2( f.v[1], f.v[0] );
  f.gam      = 1.0 / sqrt( f.vmag*f.vmag + 1.0 );
}


#endif
