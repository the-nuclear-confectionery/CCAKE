#include "system_state.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <sstream>

//using namespace std;
using namespace std::cout;
using namespace std::endl;
using namespace std::string;

#include "constants.h"
#include "vector.h"
#include "tables.h"
#include "particle.h"
#include "runge_kutta.h"
#include "eos.h"
#include "io.h"

using namespace constants;


////////////////////////////////////////////////////////////////////////////////
void SystemState::initialize()  // formerly called "manualenter"
{
  double h, factor;
  double it0;
  int start, end;

  int df;

  linklist.setv(    fvisc );
  linklist.eost   = eostype;
  linklist.cevent = 0;
  std::cout << fvisc << " hydro, h=" << h <<  " dimensions=" << D
            << " dt=" << ics.dt << " QM fluc:  " << linklist.qmf << "\n";

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
      string quantityFile   = ifolder + std::string("quantityFile.h5");
      string derivativeFile = ifolder + std::string("derivFile.h5");
      std::cout << "Using BSQ Equation of State table from: "
                << quantityFile << " and " << derivativeFile << "\n";

      EOS0.init( quantityFile, derivativeFile );
    }
    else
    {
      string quantityFilename   = "EoS_Taylor_AllMu_T0_1200.dat";
      string derivativeFilename = "EoS_Taylor_AllMu_Derivatives_T0_1200.dat";
      string quantityFile       = ifolder + quantityFilename;
      string derivativeFile     = ifolder + derivativeFilename;
      std::cout << "Using BSQ Equation of State table from: "
                << quantityFile << " and " << derivativeFile << "\n";

      EOS0.init( quantityFile, derivativeFile );
    }

    EOS0.eosin( eostype );			// does nothing!
    const double freeze_out_T_at_mu_eq_0
                  = 0.15/hbarc_GeVfm;	//1/fm
    efcheck       = EOS0.efreeze( freeze_out_T_at_mu_eq_0 );
    sfcheck       = EOS0.sfreeze( freeze_out_T_at_mu_eq_0 );
    //efcheck = 0.266112/0.1973;
    //sfcheck = 2.05743;

    std::cout << "efcheck = " << efcheck*hbarc_GeVfm << " GeV/fm^3\n";
    std::cout << "sfcheck = " << sfcheck << " 1/fm^3\n";
  }
  else
  {
    std::cerr << "This EoS model not currently supported!" << std::endl;
  }

  linklist.efcheck = efcheck;
  linklist.sfcheck = sfcheck;
  linklist.fcount  = 0;
  linklist.average = 0;
  //       Start reading ICs          //

  int numpart, _Ntable3;

  //  cout << "setting up SPH" << endl;

  cout << "Initial conditions type: " << ictype << endl;

  linklist.gtyp=0;
  if ( ictype == iccing )
  {

    int count           = 1;
    linklist.ebe_folder = outf;
    vector<string>        filelist( count );

    int j               = 0;
    filelist[j]         = ic + "/ic0.dat"; // only doing single event
    linklist.filenames  = filelist;
    linklist.fcount     = count;
    linklist.fnum       = linklist.start;

    // already done
    //readICs_iccing(linklist.filenames[0], _Ntable3, _p, factor, efcheck, numpart, EOS0);

    ////////////////////////////////////////////////////////////////////////////
    // INITIALIZE PARTICLES
    Particle::set_equation_of_state( EOS0 );

    // initialize 0th particle
    //_p[0].start(eostype, EOS0);

    // assume initial conditions have been read in from file
    

    linklist.initialize( it0, _Ntable3, h, particles, ics.dt, numpart );

    cout << "number of sph particles=" << _Ntable3 << endl;
    linklist.gtyp=6;

  }

  /*if ( ictype == iccing )
  {
    linklist.updateIC();
    cout << "bsq optimization done" << endl;
    linklist.bsqsvfreezeset();
  }*/

  // formerly bsqsv_set in this loop
  for (auto & p : particles)
  {
    double gg = p.gamcalc();
    p.g2      = gg*gg;
    p.shv33   = 0.0;
  }

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

  return;
}


void SystemState::check_BSQ_energy_conservation()
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

////////////////////////////////////////////////////////////////////////////////
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



//if we include the SPH over rhoB, rhoS, rhoQ
void SystemState::smooth_fields(int a, bool init_mode /*== false*/)
{
  const auto & pa    = _p[a];
  pa.sigma           = 0.0;
  pa.eta             = 0.0;
  pa.rhoB_sub        = 0.0;
  pa.rhoS_sub        = 0.0;
  pa.rhoQ_sub        = 0.0;
  int neighbor_count = 0;

  Vector<int,2> i;
  for ( i.x[0] =- 2; i.x[0] <= 2; i.x[0]++ )
  for ( i.x[1] =- 2; i.x[1] <= 2; i.x[1]++ )
  {
    int b = lead[ triToSum( dael[a] + i, size ) ];
    while ( b != -1 )
    {
      const auto & pb = _p[b];
      double kern     = kernel( pa.r - pb.r );
      pa.sigma       += pb.sigmaweight*kern;
      pa.eta         += pb.sigmaweight*pb.eta_sigma*kern;
      pa.rhoB_sub    += pb.rho_weight*pb.rhoB_an*kern;    //confirm with Jaki
      pa.rhoS_sub    += pb.rho_weight*pb.rhoS_an*kern;    //confirm with Jaki
      pa.rhoQ_sub    += pb.rho_weight*pb.rhoQ_an*kern;    //confirm with Jaki

      if (kern>0.0) neighbor_count++;
      if (false)
        std::cout << __FUNCTION__ << "(SPH particle == " << a << " ): "
        << b << "   " << pa.r
        << "   " << pa.sigma
        << "   " << pa.eta
        << "   " << pb.r
        << "   " << pb.sigmaweight
        << "   " << pb.eta_sigma
        << "   " << pb.rhoB_an
        << "   " << pa.rhoB_sub
        << "   " << pb.rhoS_an
        << "   " << pa.rhoS_sub
        << "   " << pb.rhoQ_an
        << "   " << pa.rhoQ_sub
        << "   " << kern << std::endl;

      b = link[b];
    }
  }

  //cout << "Check neighbor count: " << a << "   " << neighbor_count << endl;

  // reset total B, S, and Q charge of each SPH particle to reflect
  // smoothing from kernel function (ONLY ON FIRST TIME STEP)
  //cout << "-----------------------------------------------------------------" << endl;
  if ( init_mode )
  {
    //cout << "BEFORE: " << a << "   " << pa.B << "   "
    //    << pa.S << "   " << pa.Q << endl;
    //cout << pa.rho_weight << "   " << pa.rhoB_an << "   "
    //    << pa.rhoS_an << "   " << pa.rhoQ_an << endl;
    pa.B = pa.rhoB_sub * pa.rho_weight;
    pa.S = pa.rhoS_sub * pa.rho_weight;
    pa.Q = pa.rhoQ_sub * pa.rho_weight;
    //cout << "AFTER: " << a << "   " << pa.B << "   "
    //    << pa.S << "   " << pa.Q << endl;
    //cout << pa.rho_weight << "   " << pa.rhoB_sub << "   "
    //    << pa.rhoS_sub << "   " << pa.rhoQ_sub << endl;
    //cout << "-----------------------------------------------------------------" << endl;
  }

  return;
}




void SystemState::smooth_gradients( int a, double tin, int & count )
{
  const auto & pa    = _p[a];

  pa.gradP     = 0.0;
  pa.gradBulk  = 0.0;
//  pa.gradrhoB  = 0.0;
//  pa.gradrhoS  = 0.0;
//  pa.gradrhoQ  = 0.0;
  pa.gradV     = 0.0;
  pa.gradshear = 0.0;
  pa.divshear  = 0.0;

  Vector<int,2> i;

  if ( pa.btrack != -1 ) pa.btrack = 0;

  double rdis = 0;

  for ( i.x[0] =- 2; i.x[0] <= 2; i.x[0]++ )
  for ( i.x[1] =- 2; i.x[1] <= 2; i.x[1]++ )
  {

    int b=lead[ triToSum( dael[a] + i, size ) ];

    while( b != -1 )
    {
      const auto & pb          = _p[b];

      Vector<double,2> gradK   = gradKernel( pa.r - pb.r,
                                  static_cast<bool>( a == 30 && b == 43 ) );
      Vector<double,2> va      = rowp1(0, pa.shv);
      Vector<double,2> vb      = rowp1(0, pb.shv);
      Matrix<double,2,2> vminia, vminib;

      mini(vminia, pa.shv);
      mini(vminib, pb.shv);

      double sigsqra           = 1.0/(pa.sigma*pa.sigma);
      double sigsqrb           = 1.0/(pb.sigma*pb.sigma);
      Vector<double,2> sigsigK = pb.sigmaweight * pa.sigma * gradK;

      pa.gradP                += ( sigsqrb*pb.eos.p()
                                  + sigsqra*pa.eos.p() ) * sigsigK;

      if ( ( ( Norm( pa.r - pb.r ) / _h ) <= 2 ) && ( a != b ) )
      {
        if ( pa.btrack != -1 ) pa.btrack++;
        if ( pa.btrack ==  1 ) rdis = Norm(pa.r-pb.r)/_h;
      }

      pa.gradBulk             += ( pb.Bulk/pb.sigma/pb.gamma
                                    + pa.Bulk/pa.sigma/pa.gamma)/tin*sigsigK;
      //pa.gradrhoB             += ( pb.rhoB/pb.sigma/pb.gamma
      //                            + pa.rhoB/pa.sigma/pa.gamma)/tin*sigsigK;
      //pa.gradrhoS             += ( pb.rhoS/pb.sigma/pb.gamma
      //                            + pa.rhoS/pa.sigma/pa.gamma)/tin*sigsigK;
      //pa.gradrhoQ             += ( pb.rhoQ/pb.sigma/pb.gamma
      //                            + pa.rhoQ/pa.sigma/pa.gamma)/tin*sigsigK;
      pa.gradV                += pb.sigmaweight*( pb.v -  pa.v )*gradK/pa.sigma;

      pa.gradshear            += inner(sigsigK, pa.v)*( sigsqrb*vb + sigsqra*va );
      pa.divshear             += sigsqrb*sigsigK*transpose(vminib)
                                  + sigsqra*sigsigK*transpose(vminia);

      if ( isnan( pa.gradP.x[0] ) )
      {
        cout << "gradP stopped working" << endl;
        cout << t <<" "  << pa.gradP << " " << a << " " << b << endl;
        cout << pb.sigmaweight << " " << pa.sigma << " " << pb.eos.p() << endl;
        cout << Size << " " << pb.eos.s() << " " << pa.eos.s() << endl;

        cout << pa.r << endl;
        cout << pb.r << endl;
        cout << kernel( pa.r - pb.r ) << endl;
      }
      else if ( isnan( pa.gradP.x[1] ) )
        cout << "1 " << gradPressure_weight(a, b)
             << " " << a << " " << b << endl;
      else if ( isnan( pa.gradP.x[2] ) )
        cout << "2 " << gradPressure_weight(a, b)
             << " " << a << " " << b << endl;

      b=link[b];
    }
  }

  if ( ( pa.btrack == 1 )
        && ( ( pa.eos.T()*197.3 ) >= 150 ) )
    pa.frz2.t=tin;
  else if ( ( pa.btrack == 0 )
            && ( ( pa.eos.T()*197.3 ) >= 150 )
            && ( pa.Freeze < 4 ) )
    cout << "Missed " << a << " " << tin << "  "
         << pa.eos.T()*197.3 << " "
         << rdis << " " << cfon <<  endl;

  return;
}




void SystemState::bsqsvconservation()
{
    bsqsvconservation_E();
    Etot  = E + Ez;
    Eloss = (E0-Etot)/E0*100;
    rk2   = 0;
}

void SystemState::conservation_entropy()
{
  S=0.0;

  for (int i=0; i<_n; i++)
  {
    S += particles[i].eta_sigma*particles[i].sigmaweight;
    if (i==0)
    std::cout << "\t\t --> " << i << "   " << particles[i].eta_sigma << "   "
              << particles[i].sigmaweight << "   " << S << endl;
  }

  if (first==1)
    S0=S;
}


// =============================================
// function to check conservation of B, S, and Q
void SystemState::conservation_BSQ()
{
    Btotal = 0.0;
    Stotal = 0.0;
    Qtotal = 0.0;

    for (int i=0; i<_n; i++)
	{
        //Btotal += particles[i].B;
        //Stotal += particles[i].S;
        //Qtotal += particles[i].Q;
        Btotal += particles[i].rhoB_sub*particles[i].rho_weight;
        Stotal += particles[i].rhoS_sub*particles[i].rho_weight;
        Qtotal += particles[i].rhoQ_sub*particles[i].rho_weight;
    }

    if (first==1)
    {
        Btotal0 = Btotal;
        Stotal0 = Stotal;
        Qtotal0 = Qtotal;
    }
	return;
}




void SystemState::bsqsvconservation_E()
{

    E=0.;
    for (int i=0; i<_n; i++)
    {
      const auto & p = particles[i];

        E += ( p.C*p.g2 - p.eos.p() - p.bigPI + p.shv.x[0][0] )
              / p.sigma*p.sigmaweight*t;
        if (i==0)
          std::cout << "E: " << i << "   " << t
              << "   " << p.eos.T()
              << "   " << p.EOSe()
              << "   " << p.C
              << "   " << p.g2
              << "   " << p.eos.p()
              << "   " << p.bigPI
              << "   " << p.shv.x[0][0]
              << "   " << p.sigma
              << "   " << p.sigmaweight << endl;    }

    if (first==1)
    {
      first=0;
      E0=E;
    }
}




void SystemState::bsqsvconservation_Ez()
{
  dEz=0.;

  double t2=t*t;
  for (int i=0; i<_n; i++)
  {
    const auto & p = particles[i];

    dEz += ( p.eos.p() + p.bigPI + p.shv33*t2 ) / p.sigma*p.sigmaweight;

    if (false)
      std::cout << "dEz: " << i << "   " << t
        << "   " << p.eos.p()
        << "   " << p.bigPI
        << "   " << p.shv33*t2
        << "   " << p.sigma
        << "   " << p.sigmaweight << endl;
  }
}


void SystemState::setshear()
{
    for ( auto & p : particles ) p.sets(t*t);
}
