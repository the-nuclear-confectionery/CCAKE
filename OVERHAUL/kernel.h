#ifndef KERNEL_H
#define KERNEL_H

#include <stdio.h>
#include <math.h>
#include <cmath>
#include <vector>

#include "vector.h"
#include "matrix.h"
#include "mathdef.h"

namespace kernel
{
  extern double knorm, knorm2, kgrad, kgrad2;

  void set_kernel_parameters( double h );
  double kernel( double q );
  double kernel( const Vector<double,2> & a, double h );
  Vector<double,2> gradKernel( const Vector<double,2> & a, double h );
  Vector<double,2> gradKernel( const Vector<double,2> & a, double r, double h );
}

#endif