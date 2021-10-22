#include "system_state.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <sstream>

using namespace std;

#include "vector.h"
#include "tables.h"
#include "particle.h"
#include "runge_kutta.h"
#include "eos.h"
#include "io.h"


////////////////////////////////////////////////////////////////////////////////
void SystemState::initialize()  // formerly called "manualenter"
{
  double h,factor;
  double it0;
  int start,end;


  int df;

  linklist.setv(fvisc);
  linklist.eost=eostype;
  linklist.cevent=0;
  cout << fvisc << " hydro, h=" << h <<  " dimensions=" << D << " dt=" << ics.dt
        << " QM fluc:  " <<  linklist.qmf << "\n";

  //////////////////////////////////////////////////////////////////////////////
  // SET EQUATION OF STATE
  // rewrite by C. Plumberg: allow for different EOS format if using BSQ
  double efcheck = 0.0, sfcheck = 0.0;
  eos EOS0;	// now declared globally
  if ( linklist.visc == 4 )	//if we're running BSQ (table is only option)
  {
    bool using_HDF = false;
    if (using_HDF)
    {
      string quantityFile = ifolder + std::string("quantityFile.h5");
      string derivativeFile = ifolder + std::string("derivFile.h5");
      std::cout << "Using BSQ Equation of State table from: "
      << quantityFile << " and " << derivativeFile << "\n";

      EOS0.init( quantityFile, derivativeFile );
    }
    else
    {
      //string quantityFile = ifolder + std::string("quantityFile.dat");
      //string derivativeFile = ifolder + std::string("derivFile.dat");
      //string quantityFile = ifolder + std::string("full_EoS_Taylor_AllMu.dat");
      //string derivativeFile = ifolder + std::string("full_EoS_Taylor_AllMu_Derivatives.dat");
      string quantityFile = ifolder + std::string("EoS_Taylor_AllMu_T0_1200.dat");
      string derivativeFile = ifolder + std::string("EoS_Taylor_AllMu_Derivatives_T0_1200.dat");
      //string quantityFile = ifolder + std::string("EoS_Taylor_AllMu_T30_800.dat");
      //string derivativeFile = ifolder + std::string("EoS_Taylor_AllMu_Derivatives_T30_800.dat");
      //string quantityFile = ifolder + std::string("EoS_Taylor_AllMu_T30_800_dense.dat");
      //string derivativeFile = ifolder + std::string("EoS_Taylor_AllMu_Derivatives_T30_800_dense.dat");
      std::cout << "Using BSQ Equation of State table from: "
      << quantityFile << " and " << derivativeFile << "\n";

      EOS0.init( quantityFile, derivativeFile );
    }
    EOS0.eosin(eostype);			// does nothing!
    const double freeze_out_T_at_mu_eq_0 = 0.15/0.1973;	//1/fm
    efcheck = EOS0.efreeze(freeze_out_T_at_mu_eq_0);
    sfcheck = EOS0.sfreeze(freeze_out_T_at_mu_eq_0);
    //efcheck = 0.266112/0.1973;
    //sfcheck = 2.05743;

    std::cout << "efcheck = " << efcheck*0.1973 << " GeV/fm^3\n";
    std::cout << "sfcheck = " << sfcheck << " 1/fm^3\n";
  }
  else
  {
    std::cerr << "This EoS model not currently supported!" << std::endl;
  }

  linklist.efcheck=efcheck;
  linklist.sfcheck=sfcheck;
  linklist.fcount=0;
  linklist.average=0;
  //       Start reading ICs          //
  //Particle<D> *_p;
  int numpart, _Ntable3;

  //  cout << "setting up SPH" << endl;

  cout << "Initial conditions type: " << ictype << endl;

  linklist.gtyp=0;
  if (ictype==iccing)
  {

    int count=1;
    linklist.ebe_folder=outf;
    string *filelist;
    filelist=new string[count];

    int j=0;
    filelist[j]= ic+"/ic0.dat"; // only doing single event
    linklist.filenames=filelist;
    linklist.fcount=count;
    linklist.fnum=linklist.start;

    // already done
    //readICs_iccing(linklist.filenames[0], _Ntable3, _p, factor, efcheck, numpart, EOS0);

    ////////////////////////////////////////////////////////////////////////////
    // INITIALIZE PARTICLES
    Particle::set_equation_of_state( EOS0 );

    // initialize 0th particle
    _p[0].start(eostype, EOS0);
    linklist.setup(it0,_Ntable3,h,_p,ics.dt,numpart);

    cout << "number of sph particles=" << _Ntable3 << endl;
    linklist.gtyp=6;

  }

  if ((ictype==iccing))
  {
    linklist.updateIC();
    cout << "bsq optimization done" << endl;
    linklist.bsqsvfreezeset();
  }

  linklist.bsqsv_set();
  return;
}


















////////////////////////////////////////////////////////////////////////////////
void SystemState::BSQSimulation( double dt, LinkList & linklist )
{
  cout << "Ready to start hydrodynamics\n";
  linklist.frzc=0;
  linklist.cf=0;

  Output<2> out(linklist);

  BBMG<2> bbmg(linklist);
  bbmg.initial(linklist);
  cout << "started bbmg" << endl;

  linklist.t=linklist.t0;

  if ( linklist.qmf == 1 || linklist.qmf == 3 )
  {
    out.bsqsveprofile(linklist);
    cout << "printed first timestep" << endl;

    linklist.conservation_entropy();
    linklist.conservation_BSQ();

    cout << "t=" << linklist.t << " S=" << linklist.S 
         << " " << linklist.Btotal << " " << linklist.Stotal
         << " " << linklist.Qtotal << endl;

    if (linklist.qmf==1) exit(0);
  }
  else if(linklist.qmf==4)
  {
    out.eccout(linklist);
    cout << "eccentricity printed" << endl;
    exit(0);
  }


  cout << "Now let's do the main evolution!" << endl;
  linklist.Ez=0;

  while ((linklist.t<linklist.tend)&&(linklist.number_part<linklist.n()))
  {
    linklist.cfon=1;


    cout << "Entering here:" << endl;

    bsqrungeKutta2<2>( dt, &BSQshear<2>, linklist );
    linklist.conservation_entropy();
    linklist.conservation_BSQ();

    cout << "t=" << linklist.t << " " <<  linklist.Eloss << " " << linklist.S
         << " " << linklist.Btotal << " " << linklist.Stotal
         << " " << linklist.Qtotal <<  endl;

    out.bsqsveprofile(linklist);


    if (linklist.cf>0) out.bsqsvFOprint(linklist);

    if (linklist.qmf==3)
    {
      double tsub=linklist.t-floor(linklist.t);
      // if you add more points to print, must also change LinkList<D>::setup and multiply steps=floor(tend-t0)+1; by the extra number of print offs / 1fm/c
      if (tsub<(0.0+dt*0.99)||(tsub>=1-+dt*0.99)) // uncomment if you want to observe energydensity profile, conservation of energy or do a Gubser check
      {
        linklist.conservation_entropy();
        cout << "t=" << linklist.t << " S=" << linklist.S << endl;  // outputs time step
        out.bsqsveprofile(linklist);   // energy density profile
        cout << "eloss= " << linklist.t << " " <<  linklist.Eloss << endl;
        out.conservation(linklist); // conservation of energy
      }
    }

  }

  cout << "BSQ-SV simulation completed!" << endl;
  linklist.endEV();

  return;
}


////////////////////////////////////////////////////////////////////////////////
// BSQ+shear+bulk Equations of motion, only set up for 2+1 at the moment
void SystemState::BSQshear( LinkList & linklist )
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

  linklist.destroy();

  return;
}


void SystemState::check_BSQ_E_conservation()
{
  E=0.0;
  for (int i=0; i<_n; i++)
    E += ( _p[i].C* _p[i].g2 - _p[i].EOSp() - _p[i].bigPI + _p[i].shv.x[0][0] )
          *_p[i].sigmaweight*t/_p[i].sigma;

  if (first == 1)
  {
    first = 0;
    E0    = E;
  }

  return;
}

// =============================================
// function to check conservation of B, S, and Q
void SystemState::check_BSQ_charge_conservation()
{
  Btotal = 0.0;
  Stotal = 0.0;
  Qtotal = 0.0;

  for (int i=0; i<_n; i++)
  {
    //Btotal += _p[i].B;
    //Stotal += _p[i].S;
    //Qtotal += _p[i].Q;
    Btotal += _p[i].rhoB_sub*_p[i].rho_weight;
    Stotal += _p[i].rhoS_sub*_p[i].rho_weight;
    Qtotal += _p[i].rhoQ_sub*_p[i].rho_weight;
  }

  if (first==1)
  {
    Btotal0 = Btotal;
    Stotal0 = Stotal;
    Qtotal0 = Qtotal;
  }

	return;
}

