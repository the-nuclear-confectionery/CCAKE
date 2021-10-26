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
#include "linklist.h"
#include "system_state.h"
#include "runge_kutta.h"
#include "eos.h"

using namespace constants;

// Constructors and destructors.
  EquationsOfMotion::EquationsOfMotion(){}
  EquationsOfMotion::~EquationsOfMotion(){}

////////////////////////////////////////////////////////////////////////////////
// The structure here is temporary until we set the mode for different terms, 
// equations are only set up for 2+1 at the moment
void EquationsOfMotion::BSQshear( LinkList & linklist )
{
  linklist.setshear();
  linklist.initiate();

  for (int i = 0; i < linklist.n(); i++)
  {
    const auto & p = linklist._p[i];

    int curfrz = 0; //added by Christopher Plumberg to get compilation
    linklist.bsqsvoptimization(i);    // NOT bsqsvoptimization2!!!
                                      // fix arguments accordingly!!!

    if ( (p.eta<0) || isnan(p.eta) )
    {
      cout << i <<  " neg entropy " <<  p.EOST()*hbarc << " " << p.eta << endl;
      p.eta = 0;
    }

  }

  cout << "Finished first loop over SPH particles" << endl;

  int curfrz = 0;
  for ( int i = 0; i < linklist.n(); i++ )
  {
    const auto & p = linklist._p[i];

    //  Computes gamma and velocity
    p.calcbsq( linklist.t ); //resets EOS!!

    /*N.B. - eventually extend to read in viscosities from table, etc.*/
    p.setvisc( linklist.etaconst, linklist.bvf, linklist.svf,
               linklist.zTc,      linklist.sTc, linklist.zwidth,
               linklist.visc );

    if (linklist.cfon==1)
      p.frzcheck( linklist.t, curfrz, linklist.n() );

  }

  cout << "Finished second loop over SPH particles" << endl;

  if (linklist.cfon==1)
  {
    linklist.number_part += curfrz;
    linklist.list.resize(curfrz);
  }

  int m=0;
  for(int i=0; i<linklist.n(); i++)
  {
    const auto & p = linklist._p[i];

    //      Computes gradients to obtain dsigma/dt
    linklist.bsqsvoptimization2( i, linklist.t, curfrz );

    p.dsigma_dt = -p.sigma * ( p.gradV.x[0][0] + p.gradV.x[1][1] );

    p.bsqsvsigset( linklist.t, i );

    if ( (p.Freeze==3) && (linklist.cfon==1) )
    {
      linklist.list[m++] = i;
      p.Freeze           = 4;
    }

  }

  if (linklist.rk2==1)
    linklist.bsqsvconservation();

  linklist.bsqsvconservation_Ez();


  //calculate matrix elements
  for ( int i=0; i<linklist.n(); i++ )
  {
    const auto & p = linklist._p[i];

    double gamt=1./p.gamma/p.stauRelax;
    double pre=p.eta_o_tau/2./p.gamma;
    double p1=gamt-4./3./p.sigma*p.dsigma_dt+1./linklist.t/3.;
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

    linklist._p[i].div_u      = (1./ p.gamma)*inner( p.u, p.du_dt)
                              - ( p.gamma/ p.sigma ) * p.dsigma_dt ;
    linklist._p[i].bigtheta   = p.div_u*linklist.t+p.gamma;

    Matrix <double,D,D> sub   = p.pimin + p.shv.x[0][0]*p.uu/p.g2 -1./p.gamma*p.piutot;

    p.inside                  = linklist.t*(
                                inner( -minshv+p.shv.x[0][0]*p.v, p.du_dt )
                                - con2(sub, p.gradU)
                                - p.gamma*linklist.t*p.shv33 );

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


  if (linklist.cfon==1) linklist.bsqsvfreezeout(curfrz);

  return;
}






void EquationsOfMotion::BSQ(&linklist)  // shear+bulk Equations of motion, only set up completely for 2+1 at the moment
{



    linklist.setshear();
    linklist.initiate();




    for(int i=0; i<linklist.n(); i++)
    {
		//cout << "Entering this loop: i = " << i << endl;
        int curfrz=0;//added by Christopher Plumberg to get compilation
		//cout << "Calling bsqsvoptimization for i = " << i << endl;
        linklist.bsqsvoptimization(i);    // NOT bsqsvoptimization2!!! fix arguments accordingly!!!
//		linklist._p[i].rhoB_sub /= _p[i].gamma*t0;
//		linklist._p[i].rhoS_sub /= _p[i].gamma*t0;
//		linklist._p[i].rhoQ_sub /= _p[i].gamma*t0;
		//cout << "Finished bsqsvoptimization for i = " << i << endl;

        if ((linklist._p[i].eta<0)||isnan(linklist._p[i].eta))
        {
            cout << i <<  " neg entropy " <<  linklist._p[i].EOST()*197.3   << " " << linklist._p[i].eta << endl;

            linklist._p[i].eta=0;

//             int a,b;
//             a=i;
//            double sigma = 0.;
//                Vector<int,D> ii;
//            for(ii.x[0]=-2; ii.x[0]<=2; ii.x[0]++)
//                {
//                for(ii.x[1]=-2; ii.x[1]<=2; ii.x[1]++)
//                {

//                     b=linklist.lead[linklist.triToSum(linklist.dael[a]+ii, linklist.size)];
//                     while(b!=-1 )
//                     {
//                     double kern=linklist.kernel(linklist._p[b].r-linklist._p[a].r);
//                        sigma = sigma + linklist._p[b].sigmaweight*kern;
//
//
//                       b=linklist.link[b];
//                       }
//                }
//            }
//            cout << "final entropy=" << sigma << " eta_sigma=" << linklist._p[a].eta_sigma << " r=" << linklist._p[a].r << endl;
//
//            exit(1);
        }

    }

	cout << "Finished first loop over SPH particles" << endl;


    int curfrz=0;
    for(int i=0; i<linklist.n(); i++)
    {
        //  Computes gamma and velocity

		//cout << "Calling calcbsq for particle i = " << i << endl;

        linklist._p[i].calcbsq(linklist.t); //resets EOS!!
		//cout << "Finished calcbsq for particle i = " << i << endl;
        /*N.B. - eventually extend to read in viscosities from table, etc.*/linklist._p[i].setvisc(linklist.etaconst,linklist.bvf,linklist.svf,linklist.zTc,linklist.sTc,linklist.zwidth,linklist.visc);
        if (linklist.cfon==1) linklist._p[i].frzcheck(linklist.t,curfrz,linklist.n());

    }

	cout << "Finished second loop over SPH particles" << endl;

    if (linklist.cfon==1)
    {
        linklist.number_part+=curfrz;
        linklist.list.resize(curfrz);
    }

    int m=0;
    for(int i=0; i<linklist.n(); i++)
    {
        //      Computes gradients to obtain dsigma/dt
		//cout << "Calling bsqsvoptimization2 for i = " << i << endl;
        linklist.bsqsvoptimization2(i,linklist.t,curfrz);
		//cout << "Finished bsqsvoptimization2 for i = " << i << endl;

        linklist._p[i].dsigma_dt = -linklist._p[i].sigma
									*( linklist._p[i].gradV.x[0][0]
									 + linklist._p[i].gradV.x[1][1] );

/*cout << "CHECK dsigma_dt: " << i << "   " << linklist.t << "   " << linklist._p[i].dsigma_dt
		<< "   " << linklist._p[i].sigma << "   " << linklist._p[i].gradV.x[0][0]
		<< "   " << linklist._p[i].gradV.x[1][1] << endl;*/

        linklist._p[i].bsqsvsigset(linklist.t,i);
        if ((linklist._p[i].Freeze==3)&&(linklist.cfon==1))
        {
            linklist.list[m]=i;
            linklist._p[i].Freeze=4;
            ++m;
        }

    }

	cout << "Finished third loop over SPH particles" << endl;


    if (linklist.rk2==1) linklist.bsqsvconservation();
    linklist.bsqsvconservation_Ez();


    //calculate matrix elements

constexpr int ic = 0;
constexpr bool printAll = false;
    for(int i=0; i<linklist.n(); i++)
    {
        double gamt=1./linklist._p[i].gamma/linklist._p[i].stauRelax;
        double pre=linklist._p[i].eta_o_tau/2./linklist._p[i].gamma;
        //p4=gamt-linklist._p[i].sigl*4./3.;
        double p1=gamt-4./3./linklist._p[i].sigma*linklist._p[i].dsigma_dt+1./linklist.t/3.;
        Vector<double,D>  minshv=rowp1(0,linklist._p[i].shv);
        //p2=linklist._p[i].setas*gamt;
        Matrix <double,D,D> partU=linklist._p[i].gradU+transpose(linklist._p[i].gradU);

if (i==ic || printAll)
cout << "CHECK misc1: " << i << "   " << linklist.t << "   " << gamt << "   " << linklist._p[i].sigma
		<< "   " << linklist._p[i].dsigma_dt << endl;

if (i==ic || printAll)
cout << "CHECK minshv: " << i << "   " << linklist.t << "   " << minshv << endl;

if (i==ic || printAll)
cout << "CHECK partU: " << i << "   " << linklist.t << "   " << partU << endl;

        // set the Mass and the Force
        Matrix <double,D,D> M=linklist._p[i].Msub(i);
        Vector<double,D> F=linklist._p[i].Btot*linklist._p[i].u
							+ linklist._p[i].gradshear
							- ( linklist._p[i].gradP + linklist._p[i].gradBulk
								+ linklist._p[i].divshear );

if (i==ic || printAll)
cout << "CHECK M: " << i << "   " << linklist.t << "   " << M << endl;



if (i==ic || printAll)
cout << "CHECK F: " << i << "   " << linklist.t << "   " << F << "   "
		<< linklist._p[i].Btot << "   " << linklist._p[i].u << "   "
		<< linklist._p[i].gradshear << "   " << linklist._p[i].gradP << "   "
		<< linklist._p[i].gradBulk << "   " << linklist._p[i].divshear << endl;

        // shear contribution
        F+=pre*linklist._p[i].v*partU+p1*minshv;

if (i==ic || printAll)
cout << "CHECK F(again): " << i << "   " << linklist.t << "   " << F << "   "
		<< pre << "   " << linklist._p[i].v << "   " << partU << "   "
		<< p1 << "   " << minshv << endl;

//if (linklist.t > 1.8) exit(8);

        double det=deter(M);


if (i==ic || printAll)
cout << "CHECK det: " << i << "   " << linklist.t << "   " << M << "   " << det << endl;


        Matrix <double,D,D> MI;
        MI.x[0][0]=M.x[1][1]/det;
        MI.x[0][1]=-M.x[0][1]/det;
        MI.x[1][0]=-M.x[1][0]/det;
        MI.x[1][1]=M.x[0][0]/det;


if (i==ic || printAll)
cout << "CHECK MI: " << i << "   " << linklist.t << "   " << MI << endl;




        linklist._p[i].du_dt.x[0]=F.x[0]*MI.x[0][0]+F.x[1]*MI.x[0][1];
        linklist._p[i].du_dt.x[1]=F.x[0]*MI.x[1][0]+F.x[1]*MI.x[1][1];




        Matrix <double,D,D> ulpi=linklist._p[i].u*colp1(0,linklist._p[i].shv);

        double vduk=inner(linklist._p[i].v,linklist._p[i].du_dt);

        Matrix <double,D,D> Ipi=-linklist._p[i].eta_o_tau/3.*(linklist._p[i].Imat+linklist._p[i].uu)+4./3.*linklist._p[i].pimin;

        linklist._p[i].div_u = (1./ linklist._p[i].gamma)*inner( linklist._p[i].u, linklist._p[i].du_dt)
								- ( linklist._p[i].gamma/ linklist._p[i].sigma)
									* linklist._p[i].dsigma_dt ;
        linklist._p[i].bigtheta=linklist._p[i].div_u*linklist.t+linklist._p[i].gamma;

if (i==ic || printAll)
cout << "CHECK div_u: " << i
		<< "   " << linklist.t
		<< "   " << linklist._p[i].div_u
		<< "   " << linklist._p[i].gamma
		<< "   " << linklist._p[i].u
		<< "   " << linklist._p[i].du_dt
		<< "   " << inner( linklist._p[i].u, linklist._p[i].du_dt)
		<< "   " << linklist._p[i].sigma << endl;
if (i==ic || printAll)
cout << "CHECK bigtheta: " << i
		<< "   " << linklist.t
		<< "   " << linklist._p[i].bigtheta
		<< "   " << linklist._p[i].gamma << endl;

        Matrix <double,D,D> sub=linklist._p[i].pimin+linklist._p[i].shv.x[0][0]/linklist._p[i].g2*linklist._p[i].uu-1./linklist._p[i].gamma*linklist._p[i].piutot;


        linklist._p[i].inside=linklist.t*(inner((-minshv+linklist._p[i].shv.x[0][0]*linklist._p[i].v),linklist._p[i].du_dt)- con2(sub,linklist._p[i].gradU)    -      linklist._p[i].gamma*linklist.t*linklist._p[i].shv33);

if (i==ic || printAll)
std::cout << "Check inside: " << i << "   "
			<< linklist.t << "   "
			<< linklist._p[i].inside << "   "
			<< minshv << ";   "
			<< linklist._p[i].shv.x[0][0]*linklist._p[i].v << ";   "
			<< linklist._p[i].du_dt << ";   "
			<< sub << "   "
			<< linklist._p[i].gradU << ";   "
			<< linklist._p[i].gamma*linklist.t*linklist._p[i].shv33 << std::endl;



        linklist._p[i].detasigma_dt =1./linklist._p[i].sigma/linklist._p[i].EOST()
										*( -linklist._p[i].bigPI*linklist._p[i].bigtheta
											+linklist._p[i].inside);
if (i==ic || printAll)
std::cout << "Check detasigma_dt: " << i << "   "
			<< linklist.t << "   "
			<< linklist._p[i].detasigma_dt << "   "
			<< linklist._p[i].sigma << "   "
			<< linklist._p[i].EOST()*197.3 << "   "
			<< linklist._p[i].bigPI << "   "
			<< linklist._p[i].bigtheta << "   "
			<< linklist._p[i].inside << std::endl;



		// N.B. - ADD EXTRA TERMS FOR BULK EQUATION
        linklist._p[i].dBulk_dt = (-linklist._p[i].zeta/linklist._p[i].sigma*linklist._p[i].bigtheta - linklist._p[i].Bulk/linklist._p[i].gamma )/linklist._p[i].tauRelax;

        Matrix <double,D,D> ududt=linklist._p[i].u*linklist._p[i].du_dt;

		// N.B. - ADD READABLE TERM NAMES
        linklist._p[i].dshv_dt= -gamt*(linklist._p[i].pimin+linklist._p[i].setas*0.5*partU)-0.5*linklist._p[i].eta_o_tau*(ududt+transpose(ududt))+linklist._p[i].dpidtsub()-vduk*(ulpi+transpose(ulpi)+(1/linklist._p[i].gamma)*Ipi)+linklist._p[i].sigl*Ipi;

        //linklist._p[i].drhoB_dt=-linklist._p[i].rhoB*linklist._p[i].sigma*linklist._p[i].bigtheta;
        //linklist._p[i].drhoS_dt=-linklist._p[i].rhoS*linklist._p[i].sigma*linklist._p[i].bigtheta;
        //linklist._p[i].drhoQ_dt=-linklist._p[i].rhoQ*linklist._p[i].sigma*linklist._p[i].bigtheta;

    }


    if (linklist.cfon==1) linklist.bsqsvfreezeout(curfrz);


    linklist.destroy();
}










#endif