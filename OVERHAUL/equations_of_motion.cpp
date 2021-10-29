#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>

#include "constants.h"
#include "mathdef.h"
#include "vector.h"
#include "equations_of_motion.h"
#include "particle.h"
#include "system_state.h"
#include "runge_kutta.h"
#include "eos.h"

using namespace constants;

// Constructors and destructors.
  EquationsOfMotion::EquationsOfMotion(){};
  EquationsOfMotion::~EquationsOfMotion(){};

////////////////////////////////////////////////////////////////////////////////
// The structure here is temporary until we set the mode for different terms 
//which will be shear, bulk, diffusion, and coupling terms, 
//current equations are only set up for 2+1d.
void EquationsOfMotion::BSQshear( SystemState & system, SPHWorkstation & ws )
{
  ws.setshear();
  system.initialize_linklist();

  for (int i = 0; i < system.n(); i++)
  {
    auto & p = system.particles[i];

    int curfrz = 0; //added by Christopher Plumberg to get compilation
    ws.smooth_fields(i);
                                      // fix arguments accordingly!!!

    if ( (p.eta<0) || isnan(p.eta) )
    {
      cout << i <<  " neg entropy " <<  p.T()*hbarc << " " << p.eta << endl;
      p.eta = 0;
    }

  }

  cout << "Finished first loop over SPH particles" << endl;

  int curfrz = 0;
  for ( int i = 0; i < system.n(); i++ )
  {
    auto & p = system.particles[i];




    //  Computes gamma and velocity
    p.calcbsq( system.t ); //resets EOS!!




    /*N.B. - eventually extend to read in viscosities from table, etc.*/
    p.setvisc( system.etaconst, system.bvf, system.svf,
               system.zTc,      system.sTc, system.zwidth,
               system.visc );

    if (system.cfon==1)
      p.frzcheck( system.t, curfrz, system.n() );

  }


  cout << "Finished second loop over SPH particles" << endl;

  if (system.cfon==1)
  {
    system.number_part += curfrz;
    system.list.resize(curfrz);
  }

  int m=0;
  for(int i=0; i<system.n(); i++)
  {
    auto & p = system.particles[i];

    //Computes gradients to obtain dsigma/dt
    ws.smooth_gradients( i, system.t, curfrz );

    p.dsigma_dt = -p.sigma * ( p.gradV.x[0][0] + p.gradV.x[1][1] );

    p.bsqsvsigset( system.t, i );

    if ( (p.Freeze==3) && (system.cfon==1) )
    {
      system.list[m++] = i;
      p.Freeze           = 4;
    }

  }

  if (system.rk2==1)
  system.bsqsvconservation();

  system.bsqsvconservation_Ez();


  //calculate matrix elements
  for ( int i=0; i<system.n(); i++ )
  {
    auto & p = system.particles[i];

    double gamt=1./p.gamma/p.stauRelax;
    double pre=p.eta_o_tau/2./p.gamma;
    double p1=gamt-4./3./p.sigma*p.dsigma_dt+1./system.t/3.;
    Vector<double,2>  minshv=rowp1(0, p.shv);
    Matrix <double,2,2> partU = p.gradU + transpose( p.gradU );

    // set the Mass and the Force
    Matrix <double,2,2> M = p.Msub(i);
    Vector<double,2> F    = p.Btot*p.u + p.gradshear
                            - ( p.gradP + p.gradBulk + p.divshear );

    // shear contribution
    F += pre*p.v*partU + p1*minshv;

    double det=deter(M);

    Matrix <double,2,2> MI;
    MI.x[0][0]=M.x[1][1]/det;
    MI.x[0][1]=-M.x[0][1]/det;
    MI.x[1][0]=-M.x[1][0]/det;
    MI.x[1][1]=M.x[0][0]/det;


    p.du_dt.x[0]=F.x[0]*MI.x[0][0]+F.x[1]*MI.x[0][1];
    p.du_dt.x[1]=F.x[0]*MI.x[1][0]+F.x[1]*MI.x[1][1];

    Matrix <double,2,2> ulpi  = p.u*colp1(0, p.shv);

    double vduk               = inner( p.v, p.du_dt);

    Matrix <double,2,2> Ipi   = -p.eta_o_tau/3. * ( p.Imat + p.uu ) + 4./3.*p.pimin;

    system.particles[i].div_u      = (1./ p.gamma)*inner( p.u, p.du_dt)
                              - ( p.gamma/ p.sigma ) * p.dsigma_dt ;
    system.particles[i].bigtheta   = p.div_u*system.t+p.gamma;

    Matrix <double,2,2> sub   = p.pimin + (p.shv.x[0][0]/p.g2)*p.uu -1./p.gamma*p.piutot;

    p.inside                  = system.t*(
                                inner( -minshv+p.shv.x[0][0]*p.v, p.du_dt )
                                - con2(sub, p.gradU)
                                - p.gamma*system.t*p.shv33 );

    p.detasigma_dt            = 1./p.sigma/p.T()*( -p.bigPI*p.bigtheta + p.inside );
    cout << "p.inside: " << p.inside << " for particle: " << i;


    // N.B. - ADD EXTRA TERMS FOR BULK EQUATION
    p.dBulk_dt = ( -p.zeta/p.sigma*p.bigtheta - p.Bulk/p.gamma )/p.tauRelax;

    Matrix <double,2,2> ududt = p.u*p.du_dt;

    // N.B. - ADD READABLE TERM NAMES
    p.dshv_dt                 = - gamt*( p.pimin + p.setas*0.5*partU )
                               - 0.5*p.eta_o_tau*( ududt + transpose(ududt) )
                               + p.dpidtsub() + p.sigl*Ipi
                               - vduk*( ulpi + transpose(ulpi) + (1/p.gamma)*Ipi );

  }


  //if (system.cfon==1) system.bsqsvfreezeout(curfrz);

  return;
}








