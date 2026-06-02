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

// template<unsigned int D>
// KOKKOS_INLINE_FUNCTION
// double SPHkernel<D>::kernelEta(double distance, double hT )
// {
//   double norm       = knorm/Kokkos::pow(hT,2);
//   double norm2      = knorm2/Kokkos::pow(hT,2);
//   double q = distance/hT;
//   if ( q >= 2.0 )
//     return 0.0;
//   if ( q >= 1.0 )
//     return norm2*(2.0-q)*(2.0-q)*(2.0-q);

//   double qq=q*q;
//   return norm*(1 - 1.5*qq + 0.75*q*qq);
// }

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

// Generic fallback: ignore hEta, use full D-dimensional distance.
template<unsigned int D>
KOKKOS_INLINE_FUNCTION
double SPHkernel<D>::kernel(const double* r1, const double* r2,
                            double hT, double hEta)
{
  return SPHkernel<D>::kernel(SPHkernel<D>::distance(r1, r2), hT);
}

// D=3 specialization: factorized W_2D(r_perp, hT) * W_1D(|Δη|, hEta).
template<>
KOKKOS_INLINE_FUNCTION
double SPHkernel<3>::kernel(const double* r1, const double* r2,
                            double hT, double hEta)
{
  double r_perp = SPHkernel<2>::distance(r1, r2);
  double deta   = SPHkernel<1>::distance(&r1[2], &r2[2]);
  return SPHkernel<2>::kernel(r_perp, hT) * SPHkernel<1>::kernel(deta, hEta);
}

// Generic fallback: ignore hEta, use full D-dimensional distance.
template<unsigned int D>
KOKKOS_INLINE_FUNCTION
void SPHkernel<D>::gradKernel(const double* rel_sep,
                              const double* pos_a, const double* pos_b,
                              double hT, double hEta, double* grad)
{
  double r = SPHkernel<D>::distance(pos_a, pos_b);
  SPHkernel<D>::gradKernel(rel_sep, r, hT, grad);
}

// D=3 specialization: factorized gradient.
template<>
KOKKOS_INLINE_FUNCTION
void SPHkernel<3>::gradKernel(const double* rel_sep,
                              const double* pos_a, const double* pos_b,
                              double hT, double hEta, double* grad)
{
  double zero2[2] = {0.0, 0.0};
  double rel_perp[2] = {rel_sep[0], rel_sep[1]};
  double r_perp = SPHkernel<2>::distance(rel_perp, zero2);
  double deta   = SPHkernel<1>::distance(&rel_sep[2], &zero2[0]);
  double gradK_2D[2], gradK_1D;
  SPHkernel<2>::gradKernel(rel_perp, r_perp, hT, gradK_2D);
  SPHkernel<1>::gradKernel(&rel_sep[2], deta, hEta, &gradK_1D);
  double W_2D = SPHkernel<2>::kernel(r_perp, hT);
  double W_1D = SPHkernel<1>::kernel(deta, hEta);
  grad[0] = gradK_2D[0] * W_1D;
  grad[1] = gradK_2D[1] * W_1D;
  grad[2] = W_2D * gradK_1D;
}

}
