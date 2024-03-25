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

/// @brief Default constructor
template<unsigned int D>
Particle<D>::Particle() { input.s = 0.0; }

/// @brief Copy constructor
/// @tparam D The dimensionality of the problem
/// @param p The particle to be copied
/// @details Only copies the input densities, the position, velocity and shear
/// tensor.
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
