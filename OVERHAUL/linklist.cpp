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
                           vector<Particle> & particles_in,
                           double dtsave, int & numpart)
{
  t0          = it0;
  _h          = h;
  _n          = ntot;
  particles   = particles_in;
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



  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////



  // check what happens with particle separates by itself?  Where in fortran code?
  //find system boundaries

  max = particles[0].r;
  min = particles[0].r;

  for ( int i = 1; i < _n; i++ )
  for ( int j = 0; j < 2;  j++ )
  {
    if ( particles[i].r.x[j] > max.x[j] ) max.x[j] = particles[i].r.x[j];
    if ( particles[i].r.x[j] < min.x[j] ) min.x[j] = particles[i].r.x[j];
  }

  //evaluate system size

  //2*range puts extra boxes on sides of grid
  double sub = 1.0/_h;
  size       = sub*(max-min)+(2.0*range+1.0)*uni;

  //Size is the volume
  Size       = 1;

  Vector<double,2> dsub;

  // finds total volume of the system
  for ( int i = 0; i < 2; i++ )
    Size *= size.x[i];


  //dael: relates every particle with its linklist cube
  // also convert particle position to an integer

  for (int j=0; j<_n; j++)
    dael[j] = sub*(particles[j].r-min) + (1.0*range)*uni;

  //lead: relates every linklist cube with one of the particles (leader) in it
  //link: links the leader particle of one cube with the others of the same cube
  // if only one particle in cube then it is the lead

  cout << "Checking this part:" << endl;
  cout << "Size = " << Size << endl;

  lead = vector<int>(Size);

  cout << "lead.size() = " << lead.size() << endl;


  for ( int j = 0; j < Size; j++ )
    lead[j] = -1;


  for ( int k = _n-1; k >= 0; k-- )
  {
    int tt   = triToSum( dael[k], size );
    link[k]  = lead[tt];
    lead[tt] = k;
  }


  return;
}

void LinkList::initiate(vector<Particle> & particles)
{
  // check what happens with particle separates by itself?  Where in fortran code?
  //find system boundaries

  max = particles[0].r;
  min = particles[0].r;

  for ( int i = 1; i < _n; i++ )
  for ( int j = 0; j < 2;  j++ )
  {
    if ( particles[i].r.x[j] > max.x[j] ) max.x[j] = particles[i].r.x[j];
    if ( particles[i].r.x[j] < min.x[j] ) min.x[j] = particles[i].r.x[j];
  }

  //evaluate system size

  //2*range puts extra boxes on sides of grid
  double sub = 1.0/_h;
  size       = sub*(max-min)+(2.0*range+1.0)*uni;

  //Size is the volume
  Size       = 1;

  Vector<double,2> dsub;

  // finds total volume of the system
  for ( int i = 0; i < 2; i++ )
    Size *= size.x[i];


  //dael: relates every particle with its linklist cube
  // also convert particle position to an integer

  for (int j=0; j<_n; j++)
    dael[j] = sub*(particles[j].r-min) + (1.0*range)*uni;

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
