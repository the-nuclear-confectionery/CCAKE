#include <cmath>
#include <cstdio>
#include <vector>

#include "../include/eos.h"
#include "../include/formatted_output.h"
#include "../include/kernel.h"
#include "../include/linklist.h"
#include "../include/mathdef.h"
#include "../include/matrix.h"
#include "../include/particle.h"
#include "../include/stopwatch.h"
#include "../include/vector.h"

LinkList::LinkList()
{
    range = 2; //number of boxes on the sides
    for (int i=0; i<2; i++) uni(i) = 1.0;
}


void LinkList::initialize( vector<Particle> * particlesPtr_in, double h_in )
{
  h            = h_in;
  particlesPtr = particlesPtr_in;
  n_particles  = particlesPtr->size();

  link         = vector<int>(n_particles);
  dael         = vector< Vector<int,2> >(n_particles);

  kernel::set_kernel_parameters( h );

  // (re)set linklist
  reset();

  return;
}

void LinkList::reset()
{

  // min and max (x,y) coordinates of all particles
  for ( auto & p : *particlesPtr )
  for ( int j = 0; j < 2;  j++ )
  {
    if ( p.r(j) > max(j) ) max(j) = p.r(j);
    if ( p.r(j) < min(j) ) min(j) = p.r(j);
  }

  //evaluate system size [in units of h]

  //2*range puts extra boxes on sides of grid
  double inv_h = 1.0/h;
  size       = inv_h*(max-min)+(2.0*range+1.0)*uni;

  //Size is the volume
  Size       = 1;

  // finds total volume of the system
  for ( int i = 0; i < 2; i++ )
    Size *= size(i);


  //============================================================================
  //dael: relates every particle with its linklist cube
  // also convert particle position to an integer

  for (int j = 0; j < n_particles; j++)
    dael[j] = inv_h*((*particlesPtr)[j].r-min) + (1.0*range)*uni;

  //============================================================================
  //lead: relates every linklist cube with one of the particles (leader) in it
  //link: links the leader particle of one cube with the others of the same cube
  // if only one particle in cube then it is the lead

  lead = vector<int>(Size);
  for ( int j = 0; j < Size; j++ )
    lead[j] = -1;

  for ( int k = n_particles - 1; k >= 0; k-- )
  {
    int tt   = indexer( dael[k], size );
    link[k]  = lead[tt];
    lead[tt] = k;
  }




  // add vector of neighbors
  Stopwatch sw;
  sw.Start();
  
  all_neighbors.clear();
  all_neighbors.resize(n_particles);
  for (int a = 0; a < n_particles; a++)
  {
    auto & pa = all_neighbors[a];
    Vector<int,2> i;
    for ( i(0) = -2; i(0) <= 2; i(0)++ )
    for ( i(1) = -2; i(1) <= 2; i(1)++ )
    {
      int b = lead[ indexer( dael[a] + i, size ) ];
      while( b != -1 )
      {
        pa.push_back( b );
        b = link[b];
      }
    }
  }
  sw.Stop();
  formatted_output::update("set vector of neighbors in "
                            + to_string(sw.printTime()) + " s.");

}

