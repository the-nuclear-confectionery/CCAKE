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
#include "linklist.h"
#include "system_state.h"
#include "sph_workstation.h"

class EquationsOfMotion
{

public:

  EquationsOfMotion();
  ~EquationsOfMotion();
  
  void BSQshear( SystemState & system, SPHWorkstation & ws );
};

#endif