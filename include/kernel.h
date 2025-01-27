#ifndef KERNEL_H
#define KERNEL_H

#include <cmath>

#include <Cabana_Core.hpp>

namespace ccake{
/// @class SPHkernel
/// @brief A class containing static functions for the SPH kernels.
/// @details This class provides methods for calculating the SPH kernel and
/// its gradient, as well as calculating the distance between two points.
/// @tparam D The dimension of the space.
template<unsigned int D>
class SPHkernel
{
private:
  static const double knorm, knorm2, kgrad, kgrad2;

public:

  KOKKOS_FUNCTION static double kernel(double distance, double hT );
  // KOKKOS_FUNCTION static double kernelEta(double distance, double hT );
  KOKKOS_FUNCTION static void gradKernel(double const* rel_dist, double r,
                                         double h, double* grad);
  KOKKOS_FUNCTION static double distance(const double* r1, const double* r2);
};
}
#endif
