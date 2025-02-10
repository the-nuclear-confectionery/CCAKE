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

#include "../include/particle.h"

//==============================================================================
// Default constructor
Particle::Particle() { input.s = 0.0; }

//==============================================================================
// Copy constructor
Particle::Particle( const Particle& p )
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
double Particle::gamcalc() { return sqrt( Norm2(hydro.u) + 1.0 ); }

//==============================================================================
void Particle::reset_pi_tensor(double tin2)
{
  hydro.gamma    = gamcalc();
  hydro.shv(2,1) = hydro.shv(1,2);
  hydro.shv(0,1) = 1./hydro.gamma*inner(hydro.u,colp1(1,hydro.shv));
  hydro.shv(0,2) = 1./hydro.gamma*inner(hydro.u,colp1(2,hydro.shv));
  hydro.shv(1,0) = hydro.shv(0,1);
  hydro.shv(2,0) = hydro.shv(0,2);

  mini(hydro.pimin,hydro.shv);
  hydro.uu       = hydro.u*hydro.u;
  hydro.shv(0,0) = 1./hydro.gamma/hydro.gamma*con(hydro.uu,hydro.pimin);
  hydro.shv33    = (hydro.shv(0,0)-hydro.shv(1,1)-hydro.shv(2,2))/tin2;
}

