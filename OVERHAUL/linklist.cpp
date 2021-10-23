#include "vector.h"
#include "matrix.h"
#include "particle.h"
#include "mathdef.h"
#include "random.h"
#include "eos.h"
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <vector>
#include "Stopwatch.h"

// this contains functions for calling EoS directly
#include <lib.h>


LinkList::LinkList()
{
    range = 2; //number of boxes on the sides
    for (int i=0; i<2; i++) uni.x[i]=1.0;
}

LinkList::~LinkList(){}


void LinkList::destroy()
{
    delete [] lead;

}

void LinkList::setup(double it0, int ntot,double h, Particle<D> *_pin,double dtsave,int & numpart)
{
    t0=it0;
    _h=h;
    _n=ntot;
    _p=_pin;
    knorm=10/7./PI/(_h*_h);
    knorm2=knorm*0.25;
    kgrad=-10/7./PI/pow(_h,3)*3/4.;	//FIX MISSING MINUS SIGN!!!!!!  CONFIRM WITH JAKI
    kgrad2=10/7./PI/pow(_h,3)/_h;
    link=new int[_n];
    dael=new Vector<int,D>[_n];
    steps=100*(floor(tend-t0)+1);

    dt=dtsave;
    number_part=numpart;
    avgetasig=0.;
}




void LinkList::endEV()
{
    delete [] link;
    delete [] dael;
}


int LinkList::triToSum(Vector<int,2>  dael, Vector<int,2>  size)
{
    return dael.x[0] + dael.x[1]*size.x[0];

}



void LinkList::bsqsvconservation()
{

//    conservation_entropy();
    bsqsvconservation_E();
    Etot=E+Ez;
    Eloss= (E0-Etot)/E0*100;
    rk2=0;
//   cout << Eloss << "% of Energy loss at time t=" << t << endl;
//    cout << (S0-S)/S0 << "% of Entropy loss at time t=" << t << endl;

}

void LinkList::conservation_entropy()
{

    S=0.;

    for (int i=0; i<_n; i++)
    {
        S+= _p[i].eta_sigma*_p[i].sigmaweight;
      if (i==0)
          cout << "\t\t --> " << i << "   " << _p[i].eta_sigma << "   "
              << _p[i].sigmaweight << "   " << S << endl;
    }

    if (first==1)
    {

        S0=S;
    }
}


// =============================================
// function to check conservation of B, S, and Q
void LinkList::conservation_BSQ()
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




void LinkList::bsqsvconservation_E()
{

    E=0.;
    for (int i=0; i<_n; i++)
    {
        E += ( _p[i].C* _p[i].g2 - _p[i].eos.p() - _p[i].bigPI + _p[i].shv.x[0][0] )
              / _p[i].sigma*_p[i].sigmaweight*t;
if (i==0)
cout << "E: " << i << "   " << t
		<< "   " << _p[i].eos.T()
		<< "   " << _p[i].EOSe()
		<< "   " << _p[i].C
		<< "   " << _p[i].g2
		<< "   " << _p[i].eos.p()
		<< "   " << _p[i].bigPI
		<< "   " << _p[i].shv.x[0][0]
		<< "   " << _p[i].sigma
		<< "   " << _p[i].sigmaweight << endl;    }

    if (first==1)
    {
        first=0;
        E0=E;
    }
}




void LinkList::etas_set()
{

    for (int i=0; i<_n; i++) {

        _p[i].eta_sigma=_p[i].sigmaweight;
        _p[i].sigmaweight=1;
    }

}



void LinkList::bsqsv_set()
{

    for (int i=0; i<_n; i++) {
        double gg=_p[i].gamcalc();
        _p[i].g2=gg*gg;
        _p[i].shv33=0;
    }

}



void LinkList::bsqsvconservation_Ez()
{

    dEz=0.;

    double t2=t*t;
    for (int i=0; i<_n; i++) {
        dEz+=( _p[i].eos.p()+_p[i].bigPI+_p[i].shv33*t2)/_p[i].sigma*_p[i].sigmaweight;
if (false)
cout << "dEz: " << i << "   " << t
		<< "   " << _p[i].eos.p()
		<< "   " << _p[i].bigPI
		<< "   " << _p[i].shv33*t2
		<< "   " << _p[i].sigma
		<< "   " << _p[i].sigmaweight << endl;
    }

}




void LinkList::setshear()
{
    double t2=t*t;
    for (int i=0; i<_n; i++) _p[i].sets(t2);
}


void LinkList::initialize()
{
// check what happens with particle separates by itself?  Where in fortran code?
    //find system boundaries

    max=min=_p[0].r;

    for(int i=1; i<_n; i++)
    {
        for(int j=0; j<2; j++) {

            if(_p[i].r.x[j]>max.x[j]) max.x[j]=_p[i].r.x[j];

            if(_p[i].r.x[j]<min.x[j]) min.x[j]=_p[i].r.x[j];

        }
    }

    //evaluate system size

    //2*range puts extra boxes on sides of grid
    double sub=1.0/_h;

    size=sub*(max-min)+(2.0*range+1.0)*uni;

    //Size is the volume
    Size=1;

    Vector<double,2> dsub;

    // finds total volume of the system
    for(int i=0; i<2; i++) Size*=size.x[i];


    //dael: relates every particle with its linklist cube
    // also convert particle position to an integer

    for(int j=0; j<_n; j++) dael[j]=sub*(_p[j].r-min)+(1.0*range)*uni;

    //lead: relates every linklist cube with one of the particles (leader) in it
    //link: links the leader particle of one cube with the others of the same cube
    // if only one particle in cube then it is the lead

    lead = new int[Size];

    for(int j=0; j<Size; j++) lead[j]=-1;


    for(int k=_n-1; k>=0; k--)
	{
        int tt=triToSum(dael[k],size); // need to understand still... seems like it needs another coordinate
        link[k]=lead[tt];
        lead[tt]=k;
    }

    return;
}





//if we include the SPH over rhoB, rhoS, rhoQ
void LinkList::smooth_fields(int a, bool init_mode /*== false*/)
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
        std::cout << "bsqsvoptimization(SPH particle == " << a << " ): "
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




void LinkList::smooth_gradients( int a, double tin, int & count )
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
