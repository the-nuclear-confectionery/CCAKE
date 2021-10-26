#ifndef EQUATIONS_OF_MOTION_H
#define EQUATIONS_OF_MOTION_H

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
#include <stdlib.h>


#include "mathdef.h"
#include "vector.h"
#include "particle.h"
#include "system.h"
#include "system_state.h"


using std::string;
using std::vector;

class EquationsOfMotion
{

public:

  EquationsOfMotion();
  ~EquationsOfMotion();

  


private:
  // the current state of the simulation
  SystemState system;


struct Input_Parameters
{
    string IC_type; // specify initial condition type
    double h; // static SPH cutoff paramter
    double dt; // time step in fm
    double t0; // initial time in fm
    string EoS_type; // specify equation of state type
    string EoS_option; // specify specifc option for EOS
    // there should an associated EoS directory with tables
    string eta; // specificy the shear viscosity type to use
    // in transport cpefficient file
    string zeta; // specificy the bulk viscosity type to use
    // in transport cpefficient file
    double Freeze_Out_Temperature;
    string Freeze_Out_Type;
};

#endif