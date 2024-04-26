#include "kernel.h"

namespace ccake
{

template class SPHkernel<1>;
template class SPHkernel<2>;
template class SPHkernel<3>;

//1D normalizations
template<> const double SPHkernel<1>::knorm  = 2./3.;  ///< Kernel normalization for \f$ q < 1 \f$
template<> const double SPHkernel<1>::knorm2 = 1./6.;  ///< Kernel normalization for \f$ 1 < q < 2 \f$
template<> const double SPHkernel<1>::kgrad  = -1./2.; ///< Kernel gradient normalization for \f$ 1 < q < 2 \f$
template<> const double SPHkernel<1>::kgrad2 = 1.;     ///< Kernel gradient normalization for \f$ q < 1 \f$

//2D normalizations
template<> const double SPHkernel<2>::knorm  = 10./(7*M_PI);    ///< Kernel normalization for \f$ q < 1 \f$
template<> const double SPHkernel<2>::knorm2 = 5./(14.*M_PI);   ///< Kernel normalization for \f$ 1 < q < 2 \f$
template<> const double SPHkernel<2>::kgrad  = -15./(14.*M_PI); ///< Kernel gradient normalization for \f$ 1 < q < 2 \f$
template<> const double SPHkernel<2>::kgrad2 = 15./(7.*M_PI);   ///< Kernel gradient normalization for \f$ q < 1 \f$

//3D normalizations
template<> const double SPHkernel<3>::knorm  = 1./M_PI;       ///< Kernel normalization for \f$ q < 1 \f$
template<> const double SPHkernel<3>::knorm2 = 1./(4*M_PI);   ///< Kernel normalization for \f$ 1 < q < 2 \f$
template<> const double SPHkernel<3>::kgrad  = -3./(4.*M_PI); ///< Kernel gradient normalization for \f$ 1 < q < 2 \f$
template<> const double SPHkernel<3>::kgrad2 = 3./(2.*M_PI);  ///< Kernel gradient normalization for \f$ q < 1 \f$

/// @brief Calculates the SPH kernel value
/// @details This function calculates the SPH kernel for a given distance and
/// smoothing length.
/// @param distance The distance between two particles.
/// @param hT The smoothing length.
/// @return The SPH kernel value.
template<unsigned int D>
KOKKOS_INLINE_FUNCTION
double SPHkernel<D>::kernel(double distance, double hT )
{
  double norm       = knorm/Kokkos::pow(hT,D);
  double norm2      = knorm2/Kokkos::pow(hT,D);
  double q = distance/hT;
  if ( q >= 2.0 )
    return 0.0;
  if ( q >= 1.0 )
    return norm2*(2.0-q)*(2.0-q)*(2.0-q);

  double qq=q*q;
  return norm*(1 - 1.5*qq + 0.75*q*qq);
}

/// @brief Calculates the Euclidean distance between two points.
/// @details This function calculates the Euclidean distance between two points
/// in D-dimensional space.
/// @param r1 The coordinates of the first point.
/// @param r2 The coordinates of the second point.
/// @return The Euclidean distance between the two points.
template<unsigned int D>
KOKKOS_INLINE_FUNCTION
double SPHkernel<D>::distance(const double* r1, const double* r2)
{
  double d=0;
  for (int idir=0; idir<D; ++idir){
    double diff = r1[idir]-r2[idir];
    d += diff*diff;
  }
  return sqrt(d);
}

/// @brief Calculates the gradient of the kernel function
/// @details This function calculates the gradient of the SPH kernel function
/// for a given relative distance and smoothing length.
/// @tparam D The dimensionality of the simulation
/// @param rel_dist An array containing the relative distance between two 
/// particles, i.e. ra-rb
/// @param r Norm of the relative distance array above
/// @param h SPH smoothing length
/// @param grad The array to store the gradient of the kernel
template<unsigned int D>
KOKKOS_FUNCTION
void SPHkernel<D>::gradKernel(double const* rel_dist, double r, double h, double* grad )
{
  double norm_grad  = kgrad/pow(h,D+1);
  double norm_grad2 = kgrad2/pow(h,D+1);
  double q = r/h;

  if ( q >= 2.0 )
  {
    for (int idir=0; idir<D; idir++)
      grad[idir] = 0.0;
    return;
  }
  if ( q >= 1.0 )
  {
    for (int idir=0; idir<D; idir++)
      grad[idir] = (norm_grad/r)*(2.0-q)*(2.0-q)*rel_dist[idir];
    return;
  }
  for (int idir=0; idir<D; idir++)
      grad[idir] = norm_grad2*( 1.5*q - 2 )*rel_dist[idir]/h;
  return;
}
}