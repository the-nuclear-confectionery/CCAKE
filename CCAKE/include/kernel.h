#ifndef KERNEL_H
#define KERNEL_H

#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <vector>

#include <Cabana_Core.hpp>

#include "vector.h"
#include "matrix.h"
#include "mathdef.h"
#include "constants.h"

namespace kernel
{
  extern double knorm, knorm2, kgrad, kgrad2;

  //void set_kernel_parameters( double hT );
  template<unsigned int D>
  KOKKOS_INLINE_FUNCTION
  double kernel(double distance, const double hT );
  template<unsigned int D>
  KOKKOS_INLINE_FUNCTION
  void gradKernel(double const* rel_dist, double r, double h, double* grad );
  //template<unsigned int D>
  //double *gradKernel( const double *a, double r, double hT );
}

#endif