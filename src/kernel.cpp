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

// ===================== non-cubic kernels (quintic + Wendland C2/C4/C6) ======================
// The cubic spline keeps its original hand-coded path (bit-identical). The other kernels share
// a uniform form W = norm_C/h^D * f(q), q=r/h, with f(q) the dimensionless shape below.
// Norms verified by \int W d^Dx = 1; f'(q) verified vs the numerical derivative of f.

/// @brief Normalization C_D so that \int W d^Dx = 1 (D = template dimension).
template<unsigned int D>
KOKKOS_INLINE_FUNCTION
double SPHkernel<D>::norm_C(KernelType kt)
{
  const double iPI = 1.0/M_PI;
  switch (kt) {
    case KernelType::QUINTIC:
      return (D==1) ? (1.0/120.0) : (D==2) ? (7.0/478.0)*iPI : (1.0/120.0)*iPI;
    case KernelType::WENDLAND_C2:
      return (D==1) ? (5.0/8.0)   : (D==2) ? (7.0/4.0)*iPI   : (21.0/16.0)*iPI;
    case KernelType::WENDLAND_C4:
      return (D==1) ? (3.0/4.0)   : (D==2) ? (9.0/4.0)*iPI   : (495.0/256.0)*iPI;
    case KernelType::WENDLAND_C6:
      return (D==1) ? (55.0/64.0) : (D==2) ? (39.0/14.0)*iPI : (1365.0/512.0)*iPI;
    default: return 0.0; // CUBIC handled separately
  }
}

/// @brief Dimensionless shape f(q). Assumes 0 <= q < support (caller guards q>=support -> 0).
/// Wendland uses the support-2h argument u=q/2 and its 1D vs 2D/3D forms.
template<unsigned int D>
KOKKOS_INLINE_FUNCTION
double SPHkernel<D>::shape_f(double q, KernelType kt)
{
  const double u = 0.5*q;
  const double v = 1.0 - u;       // (1-u) > 0 on the Wendland support
  switch (kt) {
    case KernelType::QUINTIC: {
      const double t3 = 3.0-q; double r = t3*t3*t3*t3*t3;
      if (q < 2.0) { const double t2 = 2.0-q; r -= 6.0 *t2*t2*t2*t2*t2; }
      if (q < 1.0) { const double t1 = 1.0-q; r += 15.0*t1*t1*t1*t1*t1; }
      return r;
    }
    case KernelType::WENDLAND_C2:
      return (D==1) ? v*v*v*(1.0+3.0*u)
                    : v*v*v*v*(1.0+4.0*u);
    case KernelType::WENDLAND_C4:
      return (D==1) ? v*v*v*v*v*(1.0+5.0*u+8.0*u*u)
                    : v*v*v*v*v*v*(1.0+6.0*u+(35.0/3.0)*u*u);
    case KernelType::WENDLAND_C6:
      return (D==1) ? v*v*v*v*v*v*v*(1.0+7.0*u+19.0*u*u+21.0*u*u*u)
                    : v*v*v*v*v*v*v*v*(1.0+8.0*u+25.0*u*u+32.0*u*u*u);
    default: return 0.0;
  }
}

/// @brief f'(q) = df/dq for the non-cubic kernels.
template<unsigned int D>
KOKKOS_INLINE_FUNCTION
double SPHkernel<D>::shape_dfdq(double q, KernelType kt)
{
  const double u = 0.5*q, v = 1.0 - u;
  switch (kt) {
    case KernelType::QUINTIC: {
      const double t3 = 3.0-q; double r = -5.0*t3*t3*t3*t3;
      if (q < 2.0) { const double t2 = 2.0-q; r += 30.0 *t2*t2*t2*t2; }
      if (q < 1.0) { const double t1 = 1.0-q; r += -75.0*t1*t1*t1*t1; }
      return r;
    }
    case KernelType::WENDLAND_C2:
      return (D==1) ? -3.0*q*v*v
                    : -5.0*q*v*v*v;
    case KernelType::WENDLAND_C4:
      return (D==1) ? -7.0*u*(1.0+4.0*u)*v*v*v*v
                    : -(28.0/3.0)*u*(1.0+5.0*u)*v*v*v*v*v;
    case KernelType::WENDLAND_C6:
      return (D==1) ? v*v*v*v*v*v*(-9.0*u-54.0*u*u-105.0*u*u*u)
                    : -11.0*u*(1.0+7.0*u+16.0*u*u)*v*v*v*v*v*v*v;
    default: return 0.0;
  }
}

/// @brief Calculates the SPH kernel value
/// @details This function calculates the SPH kernel for a given distance and
/// smoothing length.
/// @param distance The distance between two particles.
/// @param hT The smoothing length.
/// @return The SPH kernel value.
template<unsigned int D>
KOKKOS_INLINE_FUNCTION
double SPHkernel<D>::kernel(double distance, double hT, KernelType kt )
{
  if (kt != KernelType::CUBIC) {
    const double q = distance/hT;
    if (q >= support(kt)) return 0.0;
    return (norm_C(kt)/Kokkos::pow(hT,D)) * shape_f(q, kt);
  }
  // --- cubic spline (legacy path, kept bit-identical) ---
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

/// @brief dW/dh for the cubic spline (used by the grad-h Omega correction).
/// @details With W = (C/h^D) f(q), q = r/h, one has
/// dW/dh = (C/h^{D+1}) [ -D f(q) - q f'(q) ], evaluated piecewise.
/// At r=0 this reduces to -D W(0,h)/h, as required.
/// @param distance |r_i - r_j|.  @param h smoothing length.  @return dW/dh.
template<unsigned int D>
KOKKOS_INLINE_FUNCTION
double SPHkernel<D>::dWdh(double distance, double h, KernelType kt)
{
  double q = distance/h;
  if (kt != KernelType::CUBIC) {
    if (q >= support(kt)) return 0.0;
    // dW/dh = -(C/h^{D+1}) [ D f(q) + q f'(q) ]
    return -(norm_C(kt)/Kokkos::pow(h,D+1))
           * ( static_cast<double>(D)*shape_f(q,kt) + q*shape_dfdq(q,kt) );
  }
  // --- cubic (legacy) ---
  if ( q >= 2.0 )
    return 0.0;
  double normp  = knorm /Kokkos::pow(h,D+1);
  double norm2p = knorm2/Kokkos::pow(h,D+1);
  if ( q >= 1.0 )
  {
    double t = 2.0 - q;
    return norm2p*( -static_cast<double>(D)*t*t*t + 3.0*q*t*t );
  }
  double qq = q*q;
  return normp*( -static_cast<double>(D)*(1.0 - 1.5*qq + 0.75*q*qq)
                 + 3.0*qq - 2.25*q*qq );
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
void SPHkernel<D>::gradKernel(double const* rel_dist, double r, double h, double* grad, KernelType kt )
{
  if (kt != KernelType::CUBIC) {
    const double q = r/h;
    if (q >= support(kt) || r <= 0.0) {
      for (int idir=0; idir<D; idir++) grad[idir] = 0.0;
      return;
    }
    // grad W = (C/h^{D+1}) f'(q) * (rel_dist / r)
    const double coef = (norm_C(kt)/pow(h,D+1)) * shape_dfdq(q,kt) / r;
    for (int idir=0; idir<D; idir++) grad[idir] = coef*rel_dist[idir];
    return;
  }
  // --- cubic spline (legacy path, kept bit-identical) ---
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
                            double hT, double hEta, KernelType kt)
{
  return SPHkernel<D>::kernel(SPHkernel<D>::distance(r1, r2), hT, kt);
}

// D=3 specialization: factorized W_2D(r_perp, hT) * W_1D(|Δη|, hEta).
template<>
KOKKOS_INLINE_FUNCTION
double SPHkernel<3>::kernel(const double* r1, const double* r2,
                            double hT, double hEta, KernelType kt)
{
  double r_perp = SPHkernel<2>::distance(r1, r2);
  double deta   = SPHkernel<1>::distance(&r1[2], &r2[2]);
  return SPHkernel<2>::kernel(r_perp, hT, kt) * SPHkernel<1>::kernel(deta, hEta, kt);
}

// Generic fallback: ignore hEta, use full D-dimensional distance.
template<unsigned int D>
KOKKOS_INLINE_FUNCTION
void SPHkernel<D>::gradKernel(const double* rel_sep,
                              const double* pos_a, const double* pos_b,
                              double hT, double hEta, double* grad, KernelType kt)
{
  double r = SPHkernel<D>::distance(pos_a, pos_b);
  SPHkernel<D>::gradKernel(rel_sep, r, hT, grad, kt);
}

// D=3 specialization: factorized gradient.
template<>
KOKKOS_INLINE_FUNCTION
void SPHkernel<3>::gradKernel(const double* rel_sep,
                              const double* pos_a, const double* pos_b,
                              double hT, double hEta, double* grad, KernelType kt)
{
  double zero2[2] = {0.0, 0.0};
  double rel_perp[2] = {rel_sep[0], rel_sep[1]};
  double r_perp = SPHkernel<2>::distance(rel_perp, zero2);
  double deta   = SPHkernel<1>::distance(&rel_sep[2], &zero2[0]);
  double gradK_2D[2], gradK_1D;
  SPHkernel<2>::gradKernel(rel_perp, r_perp, hT, gradK_2D, kt);
  SPHkernel<1>::gradKernel(&rel_sep[2], deta, hEta, &gradK_1D, kt);
  double W_2D = SPHkernel<2>::kernel(r_perp, hT, kt);
  double W_1D = SPHkernel<1>::kernel(deta, hEta, kt);
  grad[0] = gradK_2D[0] * W_1D;
  grad[1] = gradK_2D[1] * W_1D;
  grad[2] = W_2D * gradK_1D;
}

// h-weighted smoothing-length derivative. Generic fallback (D=1,2 or non-factorized
// callers): the isotropic hT*dW/dhT, so single-length grad-h is unchanged.
template<unsigned int D>
KOKKOS_INLINE_FUNCTION
double SPHkernel<D>::dWdh(const double* r1, const double* r2,
                         double hT, double hEta, KernelType kt)
{
  return hT * SPHkernel<D>::dWdh(SPHkernel<D>::distance(r1, r2), hT, kt);
}

// D=3 specialization: hT*dW/dhT + hEta*dW/dhEta for W = W_2D(r_perp,hT)*W_1D(|deta|,hEta),
// i.e. lambda*(dW/dlambda) when hT=eta_h*lambda, hEta=eta_h_eta*lambda. Uses the product
// rule across the two factors:
//   dW/dhT   = dW2Ddh(r_perp,hT) * W_1D(deta,hEta)
//   dW/dhEta = W_2D(r_perp,hT)   * dW1Ddh(deta,hEta)
template<>
KOKKOS_INLINE_FUNCTION
double SPHkernel<3>::dWdh(const double* r1, const double* r2,
                         double hT, double hEta, KernelType kt)
{
  double r_perp = SPHkernel<2>::distance(r1, r2);
  double deta   = SPHkernel<1>::distance(&r1[2], &r2[2]);
  double W2  = SPHkernel<2>::kernel(r_perp, hT, kt);
  double W1  = SPHkernel<1>::kernel(deta, hEta, kt);
  double dW2 = SPHkernel<2>::dWdh(r_perp, hT, kt);
  double dW1 = SPHkernel<1>::dWdh(deta, hEta, kt);
  return hT*dW2*W1 + hEta*W2*dW1;
}

}
