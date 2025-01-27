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
#include <omp.h>
#include <random>

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
  //phimax            = 14;
  
  //===============================================
  vjet  = 1;
  // area is taken from parameter h read in from input parameters file
  //area  = PI*pow(2.*systemPtr->h,2);

  gridx = settingsPtr->stepx;
  gridy = settingsPtr->stepy;
  //cout << "Gridx and gridy are " << gridx << "," << gridy << endl << endl;

  for (int i = 0; i < phimax; i++)
  {
    //Rjetq[i]  = 0;
    //Rjetg[i]  = 0;
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

/* 
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
      
      sph_particle.line = 0.5 * kappa * exp(z*log(settingsPtr->t0)) * exp(c*log(sph_particle.rho0)) * settingsPtr->dt; // only if initial flow=0
      //jetInfo.resize(14);
      for (int j=0; j<phimax; j++) //initializes jets at each point in grid space, over 14 directions
      {
        sph_particle.phi = phi[j];
        sph_particle.pid = j;
        //sph_particle.on  = 1;
        jetInfo.push_back(sph_particle);
      }

    }
  }
}
*/

void BBMG::initial()
{
    rho0tot = 0;
    auto& p = systemPtr->particles;
    int back_to_back = 2;
    std::random_device rd; // For true randomness
    std::mt19937 gen(rd()); // Mersenne Twister generator
    std::uniform_int_distribution<> dis(0, p.size() - 1); // Distribute over the valid range
    for (int i = 0; i < 300000; ++i)
    {
        

        int random_sph_particle = dis(gen);

        
        //Density from pressure over temperature
        double rsub = p[random_sph_particle].p() / p[random_sph_particle].T();

        rho0tot += rsub;
        if (p[random_sph_particle].T() * constants::hbarc_MeVfm > Freezeout_Temp)
        {
            field sph_particle; //field of all sph particles where we take necessary line integral info
            sph_particle.r[0] = p[random_sph_particle].r(0);
            //cout << "Positions (x) of each particle in the grid is " << p.r(0) << "\n";
            sph_particle.r[1] = p[random_sph_particle].r(1);
            sph_particle.rho0 = rsub; //Density left in terms of femtometers
            //sph_particle.sph = i;
            sph_particle.T = p[random_sph_particle].T() * constants::hbarc_MeVfm;

            double kappa = get_kappa(sph_particle.T / 1000);

            sph_particle.line = 0.5 * kappa * exp(z * log(settingsPtr->t0)) * exp(c * log(sph_particle.rho0)) * settingsPtr->dt; // only if initial flow=0
            //jetInfo.resize(14);
            //for (int j = 0; j < phimax; j++) //initializes jets at each point in grid space, over 14 directions
            std::uniform_int_distribution<> phidist(0, phimax - 1);
            for (int j = 0; j < back_to_back; j++)
            {
                sph_particle.phi = phi[phidist] + j * PI;
                sph_particle.pid = phidist + phimax * j;
                jetInfo.push_back(sph_particle);
            }
        }
    }
}


double BBMG::flow(field &f) { return f.gam*(1-f.vmag*cos(f.phi-f.vang)); }

double BBMG::efluc()
{
  int random_variable = std::rand()/RAND_MAX;
  double zeta         = random_variable*(q+2.);
  return (1.+q) / pow(q+2, 1+q) * pow(q+2.-zeta, q);
}


void BBMG::propagate()
{
  double tau  = systemPtr->t;// + settingsPtr->t0;
  //int stillon = 0;

  int countyes = 0, countno = 0;
  double Rjetnorm = 0;

  #pragma omp parallel for schedule(guided)
  for (auto& jetPropagation : jetInfo)
  {
    // propagate x,y position of jet on top of sph particles
    jetPropagation.r[0] += vjet * settingsPtr->dt * cos(jetPropagation.phi); //Flow is not here to alter the direction of the jet
    jetPropagation.r[1] += vjet * settingsPtr->dt * sin(jetPropagation.phi);
    //cout << "pid checking first: " << jetPropagation.pid << endl;

    //Move interpolation to the end of this function !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    double kappa = get_kappa(jetPropagation.T / 1000); //The /1000 here is to move temps from MeV to GeV to follow Barbara's plot, same as above
    

    //if ( /*( jetPropagation.on == 1 ) &&*/  jetPropagation.T > Freezeout_Temp ) //Can remove the if statement soon, as removal takes care of the issue 
    //{
      jetPropagation.line += kappa * exp(z*log(tau)) * exp(c*log(jetPropagation.rho)) * settingsPtr->dt * flow(jetPropagation);
      countyes++;
      //cout << "Jet directions still going: " << jetPropagation.phi << endl;
      inter( jetPropagation ); //interpolation of the area around the jet
  }
        auto condition = [this](auto& jetPropagation) {
            return jetPropagation.T <= Freezeout_Temp;
        };

        // Use std::remove_if with a lambda that captures 'condition'
        auto new_end = std::remove_if(jetInfo.begin(), jetInfo.end(),
            [this, &condition](auto& jetPropagation) {
                if (condition(jetPropagation)) {
                    jetFreezeOut.push_back(jetPropagation);
                    //jetFreezeOut.push_back(std::move(jetPropagation));
                    return true; // Mark element for removal
                }
                return false; // Keep element
            });


        for(auto& frozenJets : jetFreezeOut)
        {
          //cout << "Temperature: " << frozenJets.T << endl << "Line Integral: " << frozenJets.line << endl;
          cout << frozenJets.T << " " << frozenJets.line << " " 
               << frozenJets.rho0 << " " << frozenJets.pid << endl; // add in anything else needed
        }


        // Erase the removed elements from the source vector
        jetInfo.erase(new_end, jetInfo.end());
    int tot     = jetInfo.size();
    cout << "How many jets we have: " << tot << endl;
    int totFreeze = jetFreezeOut.size();
    cout << "How many jets froze out this timestep: " << totFreeze << endl;

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
    cout << "The magnitude of the velocity returned greater than c, something is wrong. ABORTING" << endl;
    abort();
  }
  f.vang     = atan2( f.v[1], f.v[0] );
  f.gam      = 1.0 / sqrt( 1.0 - f.vmag*f.vmag );
}


#endif
