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
using std::cout;
using std::endl;
using std::string;

#include "constants.h"
#include "vector.h"
#include "particle.h"
#include "runge_kutta.h"
#include "eos.h"
#include "io.h"

using namespace constants;


////////////////////////////////////////////////////////////////////////////////
void SystemState::set_equation_of_state( EquationOfState & eos_in )
{
  //eos = eos_in;
  eosPtr = &eos_in;
}



////////////////////////////////////////////////////////////////////////////////
void SystemState::initialize()  // formerly called "manualenter"
{
  double h, factor;
  double it0;
  int start, end;

  int df;

  linklist.setv( fvisc );
  linklist.eost   = eostype;
  linklist.cevent = 0;
  std::cout << fvisc << " hydro, h=" << h <<  " dimensions=" << D
            << " dt=" << ics.dt << " QM fluc:  " << linklist.qmf << "\n";

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

    // initialize linklist
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
///////////////////////////////////////////////////////////////////////////////
//Dekra: BSQsimulation was moved from this location to BSQhydro
///////////////////////////////////////////////////////////////////////////////
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

      pa.gradP                += ( sigsqrb*pb.eosPtr->p()
                                  + sigsqra*pa.eosPtr->p() ) * sigsigK;

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
        cout << pb.sigmaweight << " " << pa.sigma << " " << pb.eosPtr->p() << endl;
        cout << Size << " " << pb.eosPtr->s() << " " << pa.eosPtr->s() << endl;

        cout << pa.r << endl;
        cout << pb.r << endl;
        cout << kernel( pa.r - pb.r ) << endl;
      }
      else if ( isnan( pa.gradP.x[1] ) )
        cout << "1 " << gradPressure_weight(pa, pb)
             << " " << a << " " << b << endl;
      else if ( isnan( pa.gradP.x[2] ) )
        cout << "2 " << gradPressure_weight(pa, pb)
             << " " << a << " " << b << endl;

      b=link[b];
    }
  }

  if ( ( pa.btrack == 1 )
        && ( ( pa.eosPtr->T()*197.3 ) >= 150 ) )
    pa.frz2.t=tin;
  else if ( ( pa.btrack == 0 )
            && ( ( pa.eosPtr->T()*197.3 ) >= 150 )
            && ( pa.Freeze < 4 ) )
    cout << "Missed " << a << " " << tin << "  "
         << pa.eosPtr->T()*197.3 << " "
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

        E += ( p.C*p.g2 - p.eosPtr->p() - p.bigPI + p.shv.x[0][0] )
              / p.sigma*p.sigmaweight*t;
        if (i==0)
          std::cout << "E: " << i << "   " << t
              << "   " << p.eosPtr->T()
              << "   " << p.EOSe()
              << "   " << p.C
              << "   " << p.g2
              << "   " << p.eosPtr->p()
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

    dEz += ( p.eosPtr->p() + p.bigPI + p.shv33*t2 ) / p.sigma*p.sigmaweight;

    if (false)
      std::cout << "dEz: " << i << "   " << t
        << "   " << p.eosPtr->p()
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
