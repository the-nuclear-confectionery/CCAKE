#ifndef KERNEL_H
#define KERNEL_H

#include <cmath>

#include <Cabana_Core.hpp>

namespace ccake{
template<unsigned int D>
class SPHkernel
{
private:
  static const double knorm, knorm2, kgrad, kgrad2;

public:

  //void set_kernel_parameters( double hT );
  KOKKOS_FUNCTION static double kernel(double distance, double hT );
  KOKKOS_FUNCTION static void gradKernel(double const* rel_dist, double r, double h, double* grad);
  KOKKOS_FUNCTION static double distance(const double* r1, const double* r2);
  //template<unsigned int D>
  //double *gradKernel( const double *a, double r, double hT );
};
}
#endif