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
  q                 = 0; // Fluctuation parameter, subject to change
  z                 = 1; // path length dependence
  a                 = 0; //Initial jet energy dependence
  c                 = (2+z-a)/3; //medium temperature dependence
  
  //===============================================
  vjet  = 1;
  // area is taken from parameter h read in from input parameters file
  area  = PI*pow(2.*systemPtr->h,2);
  rr.resize(systemPtr->n());

  for (int i = 0; i < 14; i++)
  { // Not sure what this for yet but the quantities aren't even calculated
    Rjetq[i]  = 0;
    Rjetg[i]  = 0;
    phi[i] = i*PI/7;
  }
  //Setting final energy as a start point for the integration; This is starting in GeV
  Pfg = 10;
  Pfq = 10;

  Pfg /= constants::hbarc_GeVfm; //Converting into femtometers
  Pfq /= constants::hbarc_GeVfm;

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
      
      double kappa = get_kappa(sph_particle.T / 1000);
      
      sph_particle.line = 0.5 * kappa * pow(settingsPtr->t0, z) * pow(sph_particle.rho0, c) * settingsPtr->dt; // only if initial flow=0

      for (int j=0; j<14; j++) //initializes jets at each point in grid space, over 14 directions
      {
        sph_particle.phi = phi[j];
        sph_particle.pid = j;
        cout << "if this is not changing through 14 values we have an issue" << sph_particle.pid << endl;
        sph_particle.on  = 1;
      }
      full_sph_field.push_back(sph_particle);
    }
  }
}



double BBMG::flow(field &f) { return f.gam*(1-f.vmag*cos(f.phi-f.vang)); }


double BBMG::gftLHC(double x) { return exp(6.874017911786442 - 0.26596588700706614*pow(log(x),0.6654825884427571) - 2.677705869879314*pow(log(x),0.8020502598676339) - 2.4735502984532656*pow(log(x),0.8069542250600515) - 0.36687832133337656*pow(log(x),2.070179064516989)); } 


double BBMG::qftLHC(double x) { return exp(2.9903871818687286 - 2.0117432145703114*pow(log(x),1.0384884086567516) - 1.9187151702604879*pow(log(x),1.039887584824982) - 0.15503714201000543*pow(log(x),1.0586516925018519) - 0.15384778823017106*pow(log(x),2.0829849720841573)); }


double BBMG::efluc()
{
  int random_variable = std::rand()/RAND_MAX;
  double zeta         = random_variable*(q+2.);
  return (1.+q) / pow(q+2, 1+q) * pow(q+2.-zeta, q);
}

//KKP Fragmentation function at leading order for pions from quarks
double BBMG::fragFuncPiq(double x, double y)
{
  double lambda   = 0.088;
  double mu0      = 2;
  double sbar     = log(log(y/lambda)/log(mu0/lambda));
  double N        = 0.54610 - 0.22946*pow(sbar,1) - 0.22594*pow(sbar,2) + 0.21119*pow(sbar,3);
  double alpha    = -1.46616 - 0.45404*pow(sbar,1) - 0.12684*pow(sbar,2) + 0.27646*pow(sbar,3);
  double beta     = 1.01864 + 0.95367*pow(sbar,1) - 1.09835*pow(sbar,2) + 0.74657*pow(sbar,3);
  double gamma    = -0.01877*pow(sbar,1) + 0.02949*pow(sbar,2);
  double D        = N*pow(x,alpha)*pow(1-x,beta)*(1+gamma/x);
  return D;
}

//KKP Fragmentation function at leading order for pions from gluons
double BBMG::fragFuncPig(double x, double y)
{
  double lambda   = 0.088;
  double mu0      = 2;
  double sbar     = log(log(y/lambda)/log(mu0/lambda));
  double N        = 6.04510 - 6.61523*pow(sbar,1) - 1.64978*pow(sbar,2) + 2.68223*pow(sbar,3);
  double alpha    = -0.71378 + 0.14705*pow(sbar,1) - 1.08423*pow(sbar,2) - 0.43182*pow(sbar,3);
  double beta     = 2.92133 + 1.48429*pow(sbar,1) + 1.32887*pow(sbar,2) - 1.78696*pow(sbar,3);
  double gamma    = 0.23086*pow(sbar,1) - 0.29182*pow(sbar,2);
  double D        = N*pow(x,alpha)*pow(1-x,beta)*(1+gamma/x);
  return D;
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
    full_sph_field[i].r[0] += vjet * settingsPtr->dt * cos(full_sph_field[i].phi);
    full_sph_field[i].r[1] += vjet * settingsPtr->dt * sin(full_sph_field[i].phi);


    inter( full_sph_field[i] ); //interpolation of the field
    double kappa = get_kappa(full_sph_field[i].T / 1000); //The /1000 here is to move temps from MeV to GeV to follow Barbara's plot, same as above
    
    if ( ( full_sph_field[i].on == 1 ) && ( full_sph_field[i].T > Freezeout_Temp ) )
    {
      full_sph_field[i].line += pow(tau, z) * pow(full_sph_field[i].rho, c) * settingsPtr->dt * flow(full_sph_field[i]);
      
      //cout << "This is the value of the flow factor being multiplied: " << flow(full_sph_field[i]) << endl;
      stillon++;
    }
    else //This comes in when we drop below freezeout temp, as .on should never go to 0 on its own
    {
      full_sph_field[i].on    = 0;
      //full_sph_field[i].line += 0.5 * kappa * pow(tau,z) * pow(full_sph_field[i].rho0, c) * settingsPtr->dt; /* flow(ff[i])*/
      // Commented above out as it is still adding to the line integral, after the partons should be out of the qgp; setting to 0
      full_sph_field[i].line += 0;

      //ff[i].line *= efluc();
      // Could add in fluctuations as a multiplicative factor in the next line, like the unit converter
      P0g  = (Pfg + Cg * full_sph_field[i].line) * constants::hbarc_GeVfm; //* pow(Pfg, 1-a) 
      P0q  = (Pfq + Cq * full_sph_field[i].line) * constants::hbarc_GeVfm; //* pow(Pfq, 1-a) 
      //cout << "P0g: " << P0g << " GeV, P0q: " << P0q << " GeV" << endl;

      int jj      = full_sph_field[i].pid;
      //cout << "Value jj is taking: " << jj << endl;
      Rjetg[jj]     += pow(P0g/Pfg, 1+a) * gftLHC(P0g) / gftLHC(Pfg);
      Rjetq[jj]     += pow(P0q/Pfq, 1+a) * qftLHC(P0g) / qftLHC(Pfg); 

    }
  }

  if ( stillon == 0 )
  {
    for (int j=0; j<14; j++)
    {
      //What is going on here???
      Rjetq[j] /= rho0tot;
      Rjetg[j] /= rho0tot;
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
      double kern = kernel::kernel(rdiff);
      norm         += kern;
      double gridx = settingsPtr->stepx;
      double gridy = settingsPtr->stepy;
      //f.T      += p.T()*constants::hbarc_MeVfm*0.06*0.06*kk; //Here I changed the hardcoded grid size to a read in of the default grid from settings.
      f.T      += p.T()*constants::hbarc_MeVfm*kern;
      //cout << "Interpolated temp value is " << f.T << "\n";
      f.rho    += (p.p()/p.T())*kern;
      f.v[0]   += p.hydro.v(0)*kern;
      f.v[1]   += p.hydro.v(1)*kern;
      


      //cout << dx << " " << dy << " " << p.T()*constants::hbarc_MeVfm << " " << p.hydro.v << endl;
      
    }
  }
  //f.T       *= constants::hbarc_MeVfm;
  f.T       /= norm;
  f.rho     /= norm;
  f.v[0]    /= norm;
  f.v[1]    /= norm;
  f.vmag     = sqrt( f.v[0]*f.v[0] + f.v[1]*f.v[1] );
  if (f.vmag>1)
  {
    cout << "The magnitude of the velocity returned something greater than c, something is wrong. ABORTING" << endl;
    abort();
  }
  f.vang     = atan2( f.v[1], f.v[0] );
  f.gam      = 1.0 / sqrt( 1.0 - f.vmag*f.vmag );
}


#endif
