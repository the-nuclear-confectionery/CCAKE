#ifndef KERNEL_H
#define KERNEL_H

#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <vector>

#include <Kokkos_Macros.hpp>

#include "vector.h"
#include "matrix.h"
#include "mathdef.h"
#include "constants.h"

namespace kernel
{
  extern double knorm, knorm2, kgrad, kgrad2;

  void set_kernel_parameters( double hT );
  double kernel( double q );
  double kernel( const Vector<double,2> & a, double hT );
  template<unsigned int D>
  KOKKOS_INLINE_FUNCTION
  void gradKernel(double const* rel_dist, double r, double h, double* grad );
  //template<unsigned int D>
  //double *gradKernel( const double *a, double r, double hT );
}

#endif