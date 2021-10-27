#ifndef KERNEL_H
#define KERNEL_H

#include "vector.h"
#include "matrix.h"
#include "particle.h"
#include "mathdef.h"
#include "eos.h"
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <vector>

namespace kernel
{
  double knorm, knorm2, kgrad, kgrad2;

  double kernel( Vector<double,2> a );
  Vector<double,2> gradKernel( Vector<double,2> a );
}

#endif