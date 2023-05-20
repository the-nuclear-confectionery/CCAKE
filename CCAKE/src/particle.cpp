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

#include "particle.h"

using namespace ccake;

//Template instantiations
template class Particle<1>;
template class Particle<2>;
template class Particle<3>;

//==============================================================================
// Default constructor
template<unsigned int D>
Particle<D>::Particle() { input.s = 0.0; }

//==============================================================================
// Copy constructor
template<unsigned int D>
Particle<D>::Particle( const Particle<D>& p )
{
  r             = p.r;
  input.e       = p.input.e;
  input.rhoB    = p.input.rhoB;
  input.rhoS    = p.input.rhoS;
  input.rhoQ    = p.input.rhoQ;
  hydro.u       = p.hydro.u;
  hydro.shv     = p.hydro.shv;
  hydro.shv33   = p.hydro.shv33;
  input.s       = p.input.s;
}

//==============================================================================
template<unsigned int D>
double Particle<D>::gamcalc() { return sqrt( Norm2(hydro.u) + 1.0 ); }

//==============================================================================
template<unsigned int D>
void Particle<D>::reset_pi_tensor(double tin2)
{
  hydro.gamma    = gamcalc();
  for( unsigned int i=1; i<D+1; i++ )
    hydro.shv(0,i) = 1./hydro.gamma*inner(hydro.u,colp1(i,hydro.shv));

  for( unsigned int i=0; i<D+1; i++ )
    for( unsigned int j=i+1; j<D+1; j++ )
      hydro.shv(j,i) = hydro.shv(i,j);

  mini(hydro.pimin,hydro.shv);
  hydro.uu       = hydro.u*hydro.u;
  hydro.shv(0,0) = 1./hydro.gamma/hydro.gamma*con(hydro.uu,hydro.pimin);
  hydro.shv33    = hydro.shv(0,0);
  for ( unsigned int i=1; i<D+1; i++ )
    hydro.shv33 -= hydro.shv(i,i);
  hydro.shv33   /= tin2;
}

