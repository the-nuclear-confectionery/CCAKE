#ifndef KERNEL_H
#define KERNEL_H

#include "vector.h"
#include "matrix.h"
#include "mathdef.h"
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <vector>

namespace kernel
{
  extern double knorm, knorm2, kgrad, kgrad2;

  void set_kernel_parameters( double h );
  double kernel( Vector<double,2> a, double h );
  Vector<double,2> gradKernel( Vector<double,2> a, double h );
}

#endif