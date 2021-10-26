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

#include "mathdef.h"
#include "vector.h"
#include "equations_of_motion.h"
#include "particle.h"
#include "system.h"
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
// current equations are only set up for 2+1d.
void EquationsOfMotion::BSQshear( system & system )
{
  system.setshear();
  system.initiate();

  for (int i = 0; i < system.n(); i++)
  {
    const auto & p = system._p[i];

    int curfrz = 0; //added by Christopher Plumberg to get compilation
    system.bsqsvoptimization(i);    // NOT bsqsvoptimization2!!!
                                      // fix arguments accordingly!!!

    if ( (p.eta<0) || isnan(p.eta) )
    {
      cout << i <<  " neg entropy " <<  p.EOST()*hbarc << " " << p.eta << endl;
      p.eta = 0;
    }

  }

  cout << "Finished first loop over SPH particles" << endl;

  int curfrz = 0;
  for ( int i = 0; i < system.n(); i++ )
  {
    const auto & p = system._p[i];

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
    const auto & p = system._p[i];

    //      Computes gradients to obtain dsigma/dt
    system.bsqsvoptimization2( i, system.t, curfrz );

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
    const auto & p = system._p[i];

    double gamt=1./p.gamma/p.stauRelax;
    double pre=p.eta_o_tau/2./p.gamma;
    double p1=gamt-4./3./p.sigma*p.dsigma_dt+1./system.t/3.;
    Vector<double,D>  minshv=rowp1(0, p.shv);
    Matrix <double,D,D> partU = p.gradU + transpose( p.gradU );

    // set the Mass and the Force
    Matrix <double,D,D> M = p.Msub(i);
    Vector<double,D> F    = p.Btot*p.u + p.gradshear
                            - ( p.gradP + p.gradBulk + p.divshear );

    // shear contribution
    F += pre*p.v*partU + p1*minshv;

    double det=deter(M);

    Matrix <double,D,D> MI;
    MI.x[0][0]=M.x[1][1]/det;
    MI.x[0][1]=-M.x[0][1]/det;
    MI.x[1][0]=-M.x[1][0]/det;
    MI.x[1][1]=M.x[0][0]/det;


    p.du_dt.x[0]=F.x[0]*MI.x[0][0]+F.x[1]*MI.x[0][1];
    p.du_dt.x[1]=F.x[0]*MI.x[1][0]+F.x[1]*MI.x[1][1];

    Matrix <double,D,D> ulpi  = p.u*colp1(0, p.shv);

    double vduk               = inner( p.v, p.du_dt);

    Matrix <double,D,D> Ipi   = -p.eta_o_tau/3. * ( p.Imat + p.uu ) + 4./3.*p.pimin;

    system._p[i].div_u      = (1./ p.gamma)*inner( p.u, p.du_dt)
                              - ( p.gamma/ p.sigma ) * p.dsigma_dt ;
    system._p[i].bigtheta   = p.div_u*system.t+p.gamma;

    Matrix <double,D,D> sub   = p.pimin + p.shv.x[0][0]*p.uu/p.g2 -1./p.gamma*p.piutot;

    p.inside                  = system.t*(
                                inner( -minshv+p.shv.x[0][0]*p.v, p.du_dt )
                                - con2(sub, p.gradU)
                                - p.gamma*system.t*p.shv33 );

    p.detasigma_dt            = 1./p.sigma/p.EOST()*( -p.bigPI*p.bigtheta + p.inside );


    // N.B. - ADD EXTRA TERMS FOR BULK EQUATION
    p.dBulk_dt = ( -p.zeta/p.sigma*p.bigtheta - p.Bulk/p.gamma )/p.tauRelax;

    Matrix <double,D,D> ududt = p.u*p.du_dt;

    // N.B. - ADD READABLE TERM NAMES
    p.dshv_dt                 = - gamt*( p.pimin + p.setas*0.5*partU )
                               - 0.5*p.eta_o_tau*( ududt + transpose(ududt) )
                               + p.dpidtsub() + p.sigl*Ipi
                               - vduk*( ulpi + transpose(ulpi) + (1/p.gamma)*Ipi );

  }


  if (system.cfon==1) system.bsqsvfreezeout(curfrz);

  return;
}






void EquationsOfMotion::BSQ(&system)  // shear+bulk Equations of motion, only set up completely for 2+1 at the moment
{



    system.setshear();
    system.initiate();




    for(int i=0; i<system.n(); i++)
    {
		//cout << "Entering this loop: i = " << i << endl;
        int curfrz=0;//added by Christopher Plumberg to get compilation
		//cout << "Calling bsqsvoptimization for i = " << i << endl;
        system.bsqsvoptimization(i);    // NOT bsqsvoptimization2!!! fix arguments accordingly!!!
//		system._p[i].rhoB_sub /= _p[i].gamma*t0;
//		system._p[i].rhoS_sub /= _p[i].gamma*t0;
//		system._p[i].rhoQ_sub /= _p[i].gamma*t0;
		//cout << "Finished bsqsvoptimization for i = " << i << endl;

        if ((system._p[i].eta<0)||isnan(system._p[i].eta))
        {
            cout << i <<  " neg entropy " <<  system._p[i].EOST()*197.3   << " " << system._p[i].eta << endl;

            system._p[i].eta=0;

//             int a,b;
//             a=i;
//            double sigma = 0.;
//                Vector<int,D> ii;
//            for(ii.x[0]=-2; ii.x[0]<=2; ii.x[0]++)
//                {
//                for(ii.x[1]=-2; ii.x[1]<=2; ii.x[1]++)
//                {

//                     b=system.lead[system.triToSum(system.dael[a]+ii, system.size)];
//                     while(b!=-1 )
//                     {
//                     double kern=system.kernel(system._p[b].r-system._p[a].r);
//                        sigma = sigma + system._p[b].sigmaweight*kern;
//
//
//                       b=system.link[b];
//                       }
//                }
//            }
//            cout << "final entropy=" << sigma << " eta_sigma=" << system._p[a].eta_sigma << " r=" << system._p[a].r << endl;
//
//            exit(1);
        }

    }

	cout << "Finished first loop over SPH particles" << endl;


    int curfrz=0;
    for(int i=0; i<system.n(); i++)
    {
        //  Computes gamma and velocity

		//cout << "Calling calcbsq for particle i = " << i << endl;

        system._p[i].calcbsq(system.t); //resets EOS!!
		//cout << "Finished calcbsq for particle i = " << i << endl;
        /*N.B. - eventually extend to read in viscosities from table, etc.*/system._p[i].setvisc(system.etaconst,system.bvf,system.svf,system.zTc,system.sTc,system.zwidth,system.visc);
        if (system.cfon==1) system._p[i].frzcheck(system.t,curfrz,system.n());

    }

	cout << "Finished second loop over SPH particles" << endl;

    if (system.cfon==1)
    {
        system.number_part+=curfrz;
        system.list.resize(curfrz);
    }

    int m=0;
    for(int i=0; i<system.n(); i++)
    {
        //      Computes gradients to obtain dsigma/dt
		//cout << "Calling bsqsvoptimization2 for i = " << i << endl;
        system.bsqsvoptimization2(i,system.t,curfrz);
		//cout << "Finished bsqsvoptimization2 for i = " << i << endl;

        system._p[i].dsigma_dt = -system._p[i].sigma
									*( system._p[i].gradV.x[0][0]
									 + system._p[i].gradV.x[1][1] );

/*cout << "CHECK dsigma_dt: " << i << "   " << system.t << "   " << system._p[i].dsigma_dt
		<< "   " << system._p[i].sigma << "   " << system._p[i].gradV.x[0][0]
		<< "   " << system._p[i].gradV.x[1][1] << endl;*/

        system._p[i].bsqsvsigset(system.t,i);
        if ((system._p[i].Freeze==3)&&(system.cfon==1))
        {
            system.list[m]=i;
            system._p[i].Freeze=4;
            ++m;
        }

    }

	cout << "Finished third loop over SPH particles" << endl;


    if (system.rk2==1) system.bsqsvconservation();
    system.bsqsvconservation_Ez();


    //calculate matrix elements

constexpr int ic = 0;
constexpr bool printAll = false;
    for(int i=0; i<system.n(); i++)
    {
        double gamt=1./system._p[i].gamma/system._p[i].stauRelax;
        double pre=system._p[i].eta_o_tau/2./system._p[i].gamma;
        //p4=gamt-system._p[i].sigl*4./3.;
        double p1=gamt-4./3./system._p[i].sigma*system._p[i].dsigma_dt+1./system.t/3.;
        Vector<double,D>  minshv=rowp1(0,system._p[i].shv);
        //p2=system._p[i].setas*gamt;
        Matrix <double,D,D> partU=system._p[i].gradU+transpose(system._p[i].gradU);

if (i==ic || printAll)
cout << "CHECK misc1: " << i << "   " << system.t << "   " << gamt << "   " << system._p[i].sigma
		<< "   " << system._p[i].dsigma_dt << endl;

if (i==ic || printAll)
cout << "CHECK minshv: " << i << "   " << system.t << "   " << minshv << endl;

if (i==ic || printAll)
cout << "CHECK partU: " << i << "   " << system.t << "   " << partU << endl;

        // set the Mass and the Force
        Matrix <double,D,D> M=system._p[i].Msub(i);
        Vector<double,D> F=system._p[i].Btot*system._p[i].u
							+ system._p[i].gradshear
							- ( system._p[i].gradP + system._p[i].gradBulk
								+ system._p[i].divshear );

if (i==ic || printAll)
cout << "CHECK M: " << i << "   " << system.t << "   " << M << endl;



if (i==ic || printAll)
cout << "CHECK F: " << i << "   " << system.t << "   " << F << "   "
		<< system._p[i].Btot << "   " << system._p[i].u << "   "
		<< system._p[i].gradshear << "   " << system._p[i].gradP << "   "
		<< system._p[i].gradBulk << "   " << system._p[i].divshear << endl;

        // shear contribution
        F+=pre*system._p[i].v*partU+p1*minshv;

if (i==ic || printAll)
cout << "CHECK F(again): " << i << "   " << system.t << "   " << F << "   "
		<< pre << "   " << system._p[i].v << "   " << partU << "   "
		<< p1 << "   " << minshv << endl;

//if (system.t > 1.8) exit(8);

        double det=deter(M);


if (i==ic || printAll)
cout << "CHECK det: " << i << "   " << system.t << "   " << M << "   " << det << endl;


        Matrix <double,D,D> MI;
        MI.x[0][0]=M.x[1][1]/det;
        MI.x[0][1]=-M.x[0][1]/det;
        MI.x[1][0]=-M.x[1][0]/det;
        MI.x[1][1]=M.x[0][0]/det;


if (i==ic || printAll)
cout << "CHECK MI: " << i << "   " << system.t << "   " << MI << endl;




        system._p[i].du_dt.x[0]=F.x[0]*MI.x[0][0]+F.x[1]*MI.x[0][1];
        system._p[i].du_dt.x[1]=F.x[0]*MI.x[1][0]+F.x[1]*MI.x[1][1];




        Matrix <double,D,D> ulpi=system._p[i].u*colp1(0,system._p[i].shv);

        double vduk=inner(system._p[i].v,system._p[i].du_dt);

        Matrix <double,D,D> Ipi=-system._p[i].eta_o_tau/3.*(system._p[i].Imat+system._p[i].uu)+4./3.*system._p[i].pimin;

        system._p[i].div_u = (1./ system._p[i].gamma)*inner( system._p[i].u, system._p[i].du_dt)
								- ( system._p[i].gamma/ system._p[i].sigma)
									* system._p[i].dsigma_dt ;
        system._p[i].bigtheta=system._p[i].div_u*system.t+system._p[i].gamma;

if (i==ic || printAll)
cout << "CHECK div_u: " << i
		<< "   " << system.t
		<< "   " << system._p[i].div_u
		<< "   " << system._p[i].gamma
		<< "   " << system._p[i].u
		<< "   " << system._p[i].du_dt
		<< "   " << inner( system._p[i].u, system._p[i].du_dt)
		<< "   " << system._p[i].sigma << endl;
if (i==ic || printAll)
cout << "CHECK bigtheta: " << i
		<< "   " << system.t
		<< "   " << system._p[i].bigtheta
		<< "   " << system._p[i].gamma << endl;

        Matrix <double,D,D> sub=system._p[i].pimin+system._p[i].shv.x[0][0]/system._p[i].g2*system._p[i].uu-1./system._p[i].gamma*system._p[i].piutot;


        system._p[i].inside=system.t*(inner((-minshv+system._p[i].shv.x[0][0]*system._p[i].v),system._p[i].du_dt)- con2(sub,system._p[i].gradU)    -      system._p[i].gamma*system.t*system._p[i].shv33);

if (i==ic || printAll)
std::cout << "Check inside: " << i << "   "
			<< system.t << "   "
			<< system._p[i].inside << "   "
			<< minshv << ";   "
			<< system._p[i].shv.x[0][0]*system._p[i].v << ";   "
			<< system._p[i].du_dt << ";   "
			<< sub << "   "
			<< system._p[i].gradU << ";   "
			<< system._p[i].gamma*system.t*system._p[i].shv33 << std::endl;



        system._p[i].detasigma_dt =1./system._p[i].sigma/system._p[i].EOST()
										*( -system._p[i].bigPI*system._p[i].bigtheta
											+system._p[i].inside);
if (i==ic || printAll)
std::cout << "Check detasigma_dt: " << i << "   "
			<< system.t << "   "
			<< system._p[i].detasigma_dt << "   "
			<< system._p[i].sigma << "   "
			<< system._p[i].EOST()*197.3 << "   "
			<< system._p[i].bigPI << "   "
			<< system._p[i].bigtheta << "   "
			<< system._p[i].inside << std::endl;



		// N.B. - ADD EXTRA TERMS FOR BULK EQUATION
        system._p[i].dBulk_dt = (-system._p[i].zeta/system._p[i].sigma*system._p[i].bigtheta - system._p[i].Bulk/system._p[i].gamma )/system._p[i].tauRelax;

        Matrix <double,D,D> ududt=system._p[i].u*system._p[i].du_dt;

		// N.B. - ADD READABLE TERM NAMES
        system._p[i].dshv_dt= -gamt*(system._p[i].pimin+system._p[i].setas*0.5*partU)-0.5*system._p[i].eta_o_tau*(ududt+transpose(ududt))+system._p[i].dpidtsub()-vduk*(ulpi+transpose(ulpi)+(1/system._p[i].gamma)*Ipi)+system._p[i].sigl*Ipi;

        //system._p[i].drhoB_dt=-system._p[i].rhoB*system._p[i].sigma*system._p[i].bigtheta;
        //system._p[i].drhoS_dt=-system._p[i].rhoS*system._p[i].sigma*system._p[i].bigtheta;
        //system._p[i].drhoQ_dt=-system._p[i].rhoQ*system._p[i].sigma*system._p[i].bigtheta;

    }


    if (system.cfon==1) system.bsqsvfreezeout(curfrz);


    system.destroy();
}










#endif