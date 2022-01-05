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

////////////////////////////////////////////////////////////////////////////////
void EquationsOfMotion::set_SettingsPtr(Settings * settingsPtr_in)
{
  settingsPtr = settingsPtr_in;
}

////////////////////////////////////////////////////////////////////////////////
// The structure here is temporary until we set the mode for different terms 
// which will be shear, bulk, diffusion, and coupling terms, 
// current equations are only set up for 2+1d.
void EquationsOfMotion::BSQshear( SystemState & system, SPHWorkstation & ws )
{
  // not initial call to setshear(bool is_first_timestep)
  ws.setshear(false);
  system.reset_linklist();
  /* the above two functions should be called in BSQHydro
  (although it's not clear setshear needs to exist) */

  for (int i = 0; i < system.n(); i++)
  {
    auto & p = system.particles[i];

//if ( abs(p.r.x[0]) < 0.000001 && abs(p.r.x[1]) < 0.000001 )
//  cout << "CHECK CENTER: " << system.t << "   " << i << "   " << p.T()*hbarc << "   "
//        << p.eta/p.gamma/system.t << "   " << p.s() << endl;

    ws.smooth_fields(i);

//if ( abs(p.r.x[0]) < 0.000001 && abs(p.r.x[1]) < 0.000001 )
//  cout << "CHECK CENTER: " << system.t << "   " << i << "   " << p.T()*hbarc << "   "
//        << p.eta/p.gamma/system.t << "   " << p.s() << endl;

    if ( (p.eta<0) || isnan(p.eta) )
    {
      cout << i <<  " neg entropy " <<  p.T()*hbarc << " " << p.eta << endl;
      p.eta = 0;
    }
  }
  /* the above should be called in BSQHydro via something like
  ws.smooth_all_particle_fields() */

  cout << "Finished first loop over SPH particles" << endl;

  int curfrz = 0;
  for ( int i = 0; i < system.n(); i++ )
  {
    auto & p = system.particles[i];

    //  Computes gamma and velocity
    p.calcbsq( system.t ); //resets EOS!!
    /* would be nice to remove the above from eom,
    need to think about where to put it */

    /*N.B. - eventually extend to read in viscosities from table, etc.*/
    p.setvisc( system.etaconst, system.bvf, system.svf,
               system.zTc,      system.sTc, system.zwidth,
               system.visc );
    /* the above is obsolete when including
     transport_coefficients class */

    if ( system.cfon == 1 )
      p.frzcheck( system.t, curfrz, system.n() );
  /* not sure what cfon is but I'm sure it doesn't need to be here */


  }


  cout << "Finished second loop over SPH particles" << endl;

  if ( system.cfon == 1 )
  {
    system.number_part += curfrz;
    system.list.resize(curfrz);
  }
  /* not sure what cfon is but I'm sure it doesn't need to be here */

  int m=0;
  for ( int i=0; i<system.n(); i++ )
  {
    auto & p = system.particles[i];

    //Computes gradients to obtain dsigma/dt
    ws.smooth_gradients( i, system.t, curfrz );

    p.dsigma_dt = -p.sigma * ( p.gradV.x[0][0] + p.gradV.x[1][1] );

    p.bsqsvsigset( system.t, i );

    if ( (p.Freeze==3) && (system.cfon==1) )
    {
      system.list[m++] = i;
      p.Freeze         = 4;
    }

  }
  /* the above should probably be put into something like
  ws.smooth_all_particle_gradients() */

  if (system.rk2==1)
  system.bsqsvconservation();

  system.bsqsvconservation_Ez();

  /* conservation can definitely be called in BSQHydro */

//TRAVIS: ALL OF THE ABOVE SHOULD BE SPLIT OFF INTO DIFFERENT FUNCTIONS
// AND CALLED IN BSQHYDRO E.G. ws.smooth_gradients, system.freeze_out_check, etc

  //calculate matrix elements
constexpr int ic = 0;
constexpr bool printAll = false;
  for ( int i=0; i<system.n(); i++ )
  {
    auto & p = system.particles[i];

    double gamt = 0.0, pre = 0.0, p1 = 0.0;
    if ( settingsPtr->using_shear )
    {
      gamt = 1.0/p.gamma/p.stauRelax;
      pre  = p.eta_o_tau/p.gamma;
      p1   = gamt - 4.0/3.0/p.sigma*p.dsigma_dt + 1.0/system.t/3.0;
    }

    Vector<double,2> minshv   = rowp1(0, p.shv);
    Matrix <double,2,2> partU = p.gradU + transpose( p.gradU );

//if (i==ic || printAll)
//cout << "CHECK misc1: " << i << "   " << system.t << "   " << gamt << "   " << p.sigma
//		<< "   " << p.dsigma_dt << endl;

//if (i==ic || printAll)
//cout << "CHECK minshv: " << i << "   " << system.t << "   " << minshv << endl;

//if (i==ic || printAll)
//cout << "CHECK partU: " << i << "   " << system.t << "   " << partU << endl;


    // set the Mass and the Force
    Matrix <double,2,2> M = p.Msub(i);
    Vector<double,2> F    = p.Btot*p.u + p.gradshear
                            - ( p.gradP + p.gradBulk + p.divshear );
/* Might make more sense for M and F to be members of particle? Then the
further above loop could be done in workstation and M and F could be set
at the same time... */

//if (i==ic || printAll)
//cout << "CHECK M: " << i << "   " << system.t << "   " << M << endl;



//if (i==ic || printAll)
//cout << "CHECK F: " << i << "   " << system.t << "   " << F << "   "
//		<< p.Btot << "   " << p.u << "   "
//		<< p.gradshear << "   " << p.gradP << "   "
//		<< p.gradBulk << "   " << p.divshear << endl;

    // shear contribution
    if ( settingsPtr->using_shear )
      F += pre*p.v*partU + p1*minshv;

//if (i==ic || printAll)
//cout << "CHECK F(again): " << i << "   " << system.t << "   " << F << "   "
//		<< pre << "   " << p.v << "   " << partU << "   "
//		<< p1 << "   " << minshv << endl;

    double det=deter(M);


//if (i==ic || printAll)
//cout << "CHECK det: " << i << "   " << system.t << "   " << M << "   " << det << endl;


    Matrix <double,2,2> MI;
    MI.x[0][0]=M.x[1][1]/det;
    MI.x[0][1]=-M.x[0][1]/det;
    MI.x[1][0]=-M.x[1][0]/det;
    MI.x[1][1]=M.x[0][0]/det;
  /* This notation is still a bit weird.. but also
  MI should be a member of particle as well */

//if (i==ic || printAll)
//cout << "CHECK MI: " << i << "   " << system.t << "   " << MI << endl;


    p.du_dt.x[0]=F.x[0]*MI.x[0][0]+F.x[1]*MI.x[0][1];
    p.du_dt.x[1]=F.x[0]*MI.x[1][0]+F.x[1]*MI.x[1][1];

    Matrix <double,2,2> ulpi  = p.u*colp1(0, p.shv);

    double vduk               = inner( p.v, p.du_dt );

    Matrix <double,2,2> Ipi   = -2.0*p.eta_o_tau/3.0 * ( p.Imat + p.uu ) + 4./3.*p.pimin;

    p.div_u                   = (1./ p.gamma)*inner( p.u, p.du_dt)
                                  - ( p.gamma/ p.sigma ) * p.dsigma_dt;

    p.bigtheta                = p.div_u*system.t+p.gamma;
        /* the above lines could automaticlaly be set in particle after 
        calculating the matrix elements above */

//if (i==ic || printAll)
//cout << "CHECK div_u: " << i
//		<< "   " << system.t
//		<< "   " << p.div_u
//		<< "   " << p.gamma
//		<< "   " << p.u
//		<< "   " << p.du_dt
//		<< "   " << inner( p.u, p.du_dt)
//		<< "   " << p.sigma 
//		<< "   " << p.dsigma_dt << endl;
//if (i==ic || printAll)
//cout << "CHECK bigtheta: " << i
//		<< "   " << system.t
//		<< "   " << p.bigtheta
//		<< "   " << p.gamma << endl;

    // this term occurs in Eqs. (250) and (251) of Jaki's long notes
    // translation: pi^{ij} + pi^{00} v^i v^j - pi^{i0} v^j - pi^{0j} v^i
    Matrix <double,2,2> sub   = p.pimin + (p.shv.x[0][0]/p.g2)*p.uu -1./p.gamma*p.piutot;

    // minshv = pi^{0i}                   (i   = 1,2)
    // pimin  = pi^{ij}                   (i,j = 1,2)
    // uu     = u^i u^j                   (i,j = 1,2)
    // piu    = pi^{0i} u^j               (i,j = 1,2)
    // piutot = pi^{0i} u^j + pi^{0j} u^i (i,j = 1,2)
    // gradU  = du_i/dx^j                 (i,j = 1,2)

    if ( settingsPtr->using_shear )
      p.inside                  = system.t*(
                                inner( -minshv+p.shv.x[0][0]*p.v, p.du_dt )
                                - con2(sub, p.gradU)
                                - p.gamma*system.t*p.shv33 );

    p.detasigma_dt            = 1./p.sigma/p.T()*( -p.bigPI*p.bigtheta + p.inside );

    // N.B. - ADD EXTRA TERMS FOR BULK EQUATION
    p.dBulk_dt = ( -p.zeta/p.sigma*p.bigtheta - p.Bulk/p.gamma )/p.tauRelax;

    Matrix <double,2,2> ududt = p.u*p.du_dt;

    // N.B. - ADD READABLE TERM NAMES
    if ( settingsPtr->using_shear )
      p.dshv_dt                 = - gamt*( p.pimin + p.setas*partU )
                               - p.eta_o_tau*( ududt + transpose(ududt) )
                               + p.dpidtsub() + p.sigl*Ipi
                               - vduk*( ulpi + transpose(ulpi) + (1/p.gamma)*Ipi );

/*if (abs(p.r.x[0]+2)<1e-6 && abs(p.r.x[1])<1e-6)
{
cout << "CHECK u: " << i << "   " << p.u << endl;
cout << "CHECK du_dt: " << i << "   " << p.du_dt << endl;
cout << "CHECK gamt: " << i << "   " << gamt << endl;
cout << "CHECK pimin: " << i << "   " << p.pimin << endl;
cout << "CHECK T: " << i << "   " << p.T() << endl;
cout << "CHECK s: " << i << "   " << p.s() << endl;
cout << "CHECK p: " << i << "   " << p.p() << endl;
cout << "CHECK e: " << i << "   " << p.e() << endl;
cout << "CHECK w: " << i << "   " << p.w() << endl;
cout << "CHECK cs2: " << i << "   " << p.cs2() << endl;
cout << "CHECK setas: " << i << "   " << p.setas << endl;
cout << "CHECK partU: " << i << "   " << partU << endl;
cout << "CHECK eta_o_tau: " << i << "   " << p.eta_o_tau << endl;
cout << "CHECK sigl: " << i << "   " << p.sigl << endl;
cout << "CHECK vduk: " << i << "   " << vduk << endl;
cout << "CHECK Ipi: " << i << "   " << Ipi << endl;
cout << "CHECK stauRelax: " << i << "   " << p.stauRelax << endl;

cout << "CHECK dshv_dt: " << i
		<< " = " << p.dshv_dt
		<< " = " << - gamt*p.pimin 
		<< " = " << - p.setas*0.5*partU
		<< " = " << - 0.5*p.eta_o_tau*( ududt + transpose(ududt) )
		<< " = " << p.dpidtsub()
		<< " = " << p.sigl*Ipi
		<< " = " << - vduk*( ulpi + transpose(ulpi))
		<< " = " << - vduk* (1/p.gamma)*Ipi << endl;

exit(1);
}*/


    // Here is where we evolve the B,S,Q charge densities directly (assume ideal evolution)
    //p.drhoB_dt = -p.rhoB_sub*p.bigtheta/system.t;
    //p.drhoS_dt = -p.rhoS_sub*p.bigtheta/system.t;
    //p.drhoQ_dt = -p.rhoQ_sub*p.bigtheta/system.t;

    // assume we are working with starred (lab frame) densities per unit rapidity
    // as with sigma evolution
//    p.drhoB_dt = -p.rhoB_sub * ( p.gradV.x[0][0] + p.gradV.x[1][1] );
//    p.drhoS_dt = -p.rhoS_sub * ( p.gradV.x[0][0] + p.gradV.x[1][1] );
//    p.drhoQ_dt = -p.rhoQ_sub * ( p.gradV.x[0][0] + p.gradV.x[1][1] );


  /* all of this must be replaced for things more readable and more modular.
  Current working idea: Every term can be a separate defined function, then at run time
  terms are included based on user request by adding std:functions to a vector, each
  element storing an instance of a user specified function.. still need to work out
  the finer details... */

  }



  if (system.cfon==1)
    system.bsqsvfreezeout( curfrz );


  // keep track of which particles have left EoS grid completely
  // (reset list at end of each timestep)
  system.particles_out_of_grid.clear();
  for ( int i = 0; i < system.n(); i++ )
    if ( system.particles[i].Freeze == 5 )
      system.particles_out_of_grid.push_back( i );

  std::cout << "Summary at t = " << system.t << ": "
        << system.particles_out_of_grid.size()
        << " particles have gone out of the EoS grid." << std::endl;

  /* Not sure what any of the above does but I'm certain it can be
  done somehwere else */


  return;
}








