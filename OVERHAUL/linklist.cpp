#include "vector.h"
#include "matrix.h"
#include "particle.h"
#include "mathdef.h"
#include "eos.h"
#include "kernel.h"
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <vector>
#include "Stopwatch.h"

#include "linklist.h"

// this contains functions for calling EoS directly
#include <lib.h>


LinkList::LinkList()
{
    range = 2; //number of boxes on the sides
    for (int i=0; i<2; i++) uni.x[i]=1.0;
}

//LinkList::~LinkList(){}


void LinkList::initialize( double it0, int ntot, double h,
                           vector<Particle> * particlesPtr_in,
                           double dtsave, int & numpart)
{
  t0          = it0;
  _h          = h;
  _n          = ntot;

  cout << "Check read in" << endl;
  cout << "t0 = " << t0 << endl;
  cout << "_h = " << _h << endl;
  cout << "_n = " << _n << endl;

//    cout << "Check 1: " << particles_in[0].r.x[0] << "   " << particles_in[0].r.x[1] << endl;

  //particles   = particles_in;
  particlesPtr = particlesPtr_in;

//    cout << "Check 2: " << particles[0].r.x[0] << "   " << particles[0].r.x[1] << endl;

//    knorm       = 10/7./pi/(_h*_h);
//    knorm2      = knorm*0.25;
//    kgrad       = -10/7./pi/pow(_h,3)*3/4.;	//FIX MISSING MINUS SIGN!!!!!!  CONFIRM WITH JAKI
//    kgrad2      = 10/7./pi/pow(_h,3)/_h;
  //link        = new int[_n];
  //dael        = new Vector<int,2>[_n];
  link        = vector<int>(_n);
  dael        = vector< Vector<int,2> >(_n);
  steps       = 100*(floor(tend-t0)+1);

  kernel::set_kernel_parameters( _h );

  dt          = dtsave;
  number_part = numpart;
  avgetasig   = 0.0;


  // initialize linklist
  initiate();



  return;
}

void LinkList::initiate()
{
  // check what happens with particle separates by itself?  Where in fortran code?
  //find system boundaries

  //max = (*particlesPtr)[0].r;
  //min = (*particlesPtr)[0].r;

  //for ( int i = 1; i < _n; i++ )
  for ( auto & p : *particlesPtr )
  for ( int j = 0; j < 2;  j++ )
  {
    if ( p.r.x[j] > max.x[j] ) max.x[j] = p.r.x[j];
    if ( p.r.x[j] < min.x[j] ) min.x[j] = p.r.x[j];
  }

  //evaluate system size

  //2*range puts extra boxes on sides of grid
  double inv_h = 1.0/_h;
  size       = inv_h*(max-min)+(2.0*range+1.0)*uni;

  //Size is the volume
  Size       = 1;

  // finds total volume of the system
  for ( int i = 0; i < 2; i++ )
    Size *= size.x[i];


  //dael: relates every particle with its linklist cube
  // also convert particle position to an integer

  for (int j=0; j<_n; j++)
    dael[j] = inv_h*((*particlesPtr)[j].r-min) + (1.0*range)*uni;
  //for ( auto & p : *particlesPtr )
  //  dael.push_back( inv_h*(p.r-min) + (1.0*range)*uni );

  //lead: relates every linklist cube with one of the particles (leader) in it
  //link: links the leader particle of one cube with the others of the same cube
  // if only one particle in cube then it is the lead

  lead = vector<int>(Size);

  for ( int j = 0; j < Size; j++ )
    lead[j] = -1;


  for ( int k = _n-1; k >= 0; k-- )
  {
    int tt   = triToSum( dael[k], size );
    link[k]  = lead[tt];
    lead[tt] = k;
  }

}

int LinkList::triToSum( Vector<int,2> dael, Vector<int,2> size )
{
    return dael.x[0] + dael.x[1]*size.x[0];
}
