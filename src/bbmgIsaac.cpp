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

  gridx = settingsPtr->stepx;
  gridy = settingsPtr->stepy;
  cout << "Gridx and gridy are " << gridx << "," << gridy << endl << endl;

  for (int i = 0; i < 14; i++)
  {
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
      
      sph_particle.line = 0.5 * kappa * exp(z*log(settingsPtr->t0)) * exp(c*log(sph_particle.rho0)) * settingsPtr->dt; // only if initial flow=0
      //jetInfo.resize(14);
      for (int j=0; j<14; j++) //initializes jets at each point in grid space, over 14 directions
      {
        sph_particle.phi = phi[j];
        sph_particle.pid = j;
        //sph_particle.on  = 1;
        jetInfo.push_back(sph_particle); // Attempting...
      }
      // Would putting this vector pushback inside for loop fix everything? 
      // Seems to have fixed, leaving note here for now to know where I screwed up if it's actually wrong
      //full_sph_field.push_back(sph_particle);
    }
  }
}



double BBMG::flow(field &f) { return f.gam*(1-f.vmag*cos(f.phi-f.vang)); }



// G0s and PDFs can be removed from here and will be dealt with in post processing



double BBMG::gftLHC(double x) 
{ 
  double logx = log(x);
  double loglogx = log(logx);
  return exp(6.874017911786442 
            - 0.26596588700706614*exp(0.6654825884427571*loglogx) - 2.677705869879314*exp(0.8020502598676339*loglogx) 
            - 2.4735502984532656*exp(0.8069542250600515*loglogx) - 0.36687832133337656*exp(2.070179064516989*loglogx)); 
}
//double BBMG::gftLHC(double x) { return 1;} 


double BBMG::qftLHC(double x) 
{ 
  double logx = log(x);
  double loglogx = log(logx);
  return exp(2.9903871818687286 
          - 2.0117432145703114*exp(1.0384884086567516*loglogx) - 1.9187151702604879*exp(1.039887584824982*loglogx)
          - 0.15503714201000543*exp(1.0586516925018519*loglogx) - 0.15384778823017106*exp(2.0829849720841573*loglogx));
}
//double BBMG::qftLHC(double x) { return 1;}

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
  double sbar2    = sbar * sbar;
  double sbar3    = sbar2 * sbar;
  double N        = 0.54610 - 0.22946*sbar - 0.22594*sbar2 + 0.21119*sbar3;
  double alpha    = -1.46616 - 0.45404*sbar - 0.12684*sbar2 + 0.27646*sbar3;
  double beta     = 1.01864 + 0.95367*sbar - 1.09835*sbar2 + 0.74657*sbar3;
  double gamma    = -0.01877*sbar + 0.02949*sbar2;
  double D        = N*exp(alpha*log(x))*exp(beta*log(1-x))*(1+gamma/x);
  return D;
}

//KKP Fragmentation function at leading order for pions from gluons
double BBMG::fragFuncPig(double x, double y)
{
  double lambda   = 0.088;
  double mu0      = 2;
  double sbar     = log(log(y/lambda)/log(mu0/lambda));
  double sbar2    = sbar * sbar;
  double sbar3    = sbar2 * sbar;
  double N        = 6.04510 - 6.61523*sbar - 1.64978*sbar2 + 2.68223*sbar3;
  double alpha    = -0.71378 + 0.14705*sbar - 1.08423*sbar2 - 0.43182*sbar3;
  double beta     = 2.92133 + 1.48429*sbar + 1.32887*sbar2 - 1.78696*sbar3;
  double gamma    = 0.23086*sbar - 0.29182*sbar2;
  double D        = N*exp(alpha*log(x))*exp(beta*log(1-x))*(1+gamma/x);
  return D;
}


//Delaunay triangulation if needed
/*void BBMG::initialize_Delaunay()
{

}*/

void BBMG::propagate()
{
  double tau  = systemPtr->t + settingsPtr->t0;
  //int stillon = 0;
  int tot     = jetInfo.size();
  cout << "How many jets we have: " << tot << endl;
  int countyes = 0, countno = 0;
  double P0g = 0, P0q = 0;
  double Rjetnorm = 0;
  double g0Pfg = gftLHC(Pfg);
  double g0Pfq = qftLHC(Pfq);
  // for (int i = 0; i < tot; i++)
  // Since this is in place, get rid of the if statement and figure out how to print information
  // now that we can use system state and freeze out as examples to shift jets from one to the next.
  
  // Define the condition for moving elements
    auto condition = [this](auto& jetPropagation) {
            return jetPropagation.T <= Freezeout_Temp;
        };

        // Use std::remove_if with a lambda that captures 'condition'
        auto new_end = std::remove_if(jetInfo.begin(), jetInfo.end(),
            [this, &condition](JetPropagation& jetPropagation) {
                if (condition(jetPropagation)) {
                    jetFreezeOut.push_back(std::move(jetPropagation));
                    return true; // Mark element for removal
                }
                return false; // Keep element
            });

        // Erase the removed elements from the source vector
        jetInfo.erase(new_end, jetInfo.end());
    }
  
  
  /*jetInfo.erase( std::remove_if(
      jetInfo.begin(),
      jetInfo.end(),
      [this]( auto& jetPropagation ){ return jetPropagation.T <= Freezeout_Temp; }),
      jetInfo.end() );*/

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
      jetPropagation.line += exp(z*log(tau)) * exp(c*log(jetPropagation.rho)) * settingsPtr->dt * flow(jetPropagation);
      countyes++;
      cout << "Jet directions still going: " << jetPropagation.phi << endl;
      inter( jetPropagation ); //interpolation of the area around the jet


      //cout << "This is the value of the flow factor being multiplied: " << flow(jetPropagation) << endl;
      //stillon++;
    //}
/*    else //This comes in when we drop below freezeout temp, as .on should never go to 0 on its own
    {

      jetPropagation.line += 0;
      countno++;
      //ff[i].line *= efluc();
      // Could add in fluctuations as a multiplicative factor in the next line, like the unit converter
      P0g  = (Pfg + Cg * jetPropagation.line) * constants::hbarc_GeVfm; //* pow(Pfg, 1-a) 
      P0q  = (Pfq + Cq * jetPropagation.line) * constants::hbarc_GeVfm; //* pow(Pfq, 1-a) 
      //cout << "Initial quark jet energy is: " << P0q << " GeV " << endl << "Initial gluon jet energy is: " << P0g << " GeV" << endl;

      int jj      = jetPropagation.pid;
      
      //Following Barbara's format here from Rjet_g3 or Rjet_q3 to _g11/_q11
      Rjetg[jj]     += (exp((1+a)*log(P0g/Pfg)) * gftLHC(P0g) / g0Pfg) * jetPropagation.rho0 * gridx*gridy;
      Rjetq[jj]     += (exp((1+a)*log(P0q/Pfq)) * qftLHC(P0q) / g0Pfq) * jetPropagation.rho0 * gridx*gridy; 
      
      Rjetnorm += jetPropagation.rho0 * gridx*gridy;
      //stillon = 0;
    }
  }

  for (int j=0; j<14; j++)
  {
    //There is an experimental reason to be calculating this way
    Rjetq[j] /= Rjetnorm;
    Rjetg[j] /= Rjetnorm;
    //cout << "Checking Barbara's Rjet12" << Rjetnorm << endl;
    cout << "The averaged quark jet RAA for " << j << " is: " << Rjetq[j] << endl;
    cout << "The averaged gluon jet RAA for " << j << " is: " << Rjetg[j] << endl;
  }
  cout << "Frozen out jets: " << countno << endl << "Still going jets: " << countyes << endl;
*/}



/*double BBMG::int1(double x)
{
  double p_pi_int = Pfg/x;
  double integrand1 = 1/x * gftLHC(p_pi_int) * Rjetg[j] * fragFuncPig(x, p_pi_int);
  // CONTINUE WRITING INTEGRATION FUNCTIONS
  // The returned values should be arrays for pion RAA calculation

  return integrand1;
}

double BBMG::int2(double x)
{
  double p_pi_int = Pfq/x;
  double integrand2 = 1/x * qftLHC(p_pi_int) * Rjetq[j] * fragFuncPiq(x, p_pi_int);

  return integrand2;
}*/

double BBMG::int3(double x)
{
  double p_pi_int = Pfg/x;
  double integrand3 = 1/x * gftLHC(p_pi_int) * fragFuncPig(x, p_pi_int);

  return integrand3;
}

double BBMG::int4(double x)
{
  double p_pi_int = Pfq/x;
  double integrand4 = 1/x * qftLHC(p_pi_int) * fragFuncPiq(x, p_pi_int);

  return integrand4;
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
