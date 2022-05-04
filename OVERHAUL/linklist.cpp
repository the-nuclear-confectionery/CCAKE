#include <cmath>
#include <cstdio>
#include <vector>

#include "eos.h"
#include "kernel.h"
#include "linklist.h"
#include "mathdef.h"
#include "matrix.h"
#include "particle.h"
#include "Stopwatch.h"
#include "vector.h"

LinkList::LinkList()
{
    range = 2; //number of boxes on the sides
    for (int i=0; i<2; i++) uni(i) = 1.0;
}


void LinkList::initialize( double it0, int ntot, double h,
                           vector<Particle> * particlesPtr_in,
                           double dtsave, int & numpart )
{
  t0          = it0;
  _h          = h;
  _n          = ntot;

  cout << "Check read in" << endl;
  cout << "t0 = " << t0 << endl;
  cout << "_h = " << _h << endl;
  cout << "_n = " << _n << endl;

  particlesPtr = particlesPtr_in;

  link        = vector<int>(_n);
  dael        = vector< Vector<int,2> >(_n);
  steps       = 100*(floor(tend-t0)+1);

  kernel::set_kernel_parameters( _h );

  dt          = dtsave;

  // initialize linklist
  reset();



  return;
}

void LinkList::reset()
{
  for ( auto & p : *particlesPtr )
  for ( int j = 0; j < 2;  j++ )
  {
    if ( p.r(j) > max(j) ) max(j) = p.r(j);
    if ( p.r(j) < min(j) ) min(j) = p.r(j);
  }

  //evaluate system size

  //2*range puts extra boxes on sides of grid
  double inv_h = 1.0/_h;
  size       = inv_h*(max-min)+(2.0*range+1.0)*uni;

  //Size is the volume
  Size       = 1;

  // finds total volume of the system
  for ( int i = 0; i < 2; i++ )
    Size *= size(i);


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

int LinkList::triToSum( Vector<int,2> dael_local, Vector<int,2> size_local )
{
    return dael_local(0) + dael_local(1)*size_local(0);
}
