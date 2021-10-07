#ifndef PARTICLE_H
#define PARTICLE_H

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

struct SPH_thermo
{
  double T, muB, muS, muQ;
  double e, p, s, B, S, Q;
  double dwds, dwdB, dwdS, dwdQ;
}

class particle
{
  // Constructors and destructors.
  // particle();
  // ~particle();


  private:
    vector<double> position;  // (x,y) in fm

    SPH_thermo thermodynamic_state;


  public:



}

#endif