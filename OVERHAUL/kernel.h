#include "vector.h"
#include "matrix.h"
#include "particle.h"
#include "mathdef.h"
#include "random.h"
#include "eos.h"
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <vector>

namespace kernel
{
  double kernel( Vector<double,2> a );
  Vector<double,2> gradKernel( Vector<double,2> a );
}