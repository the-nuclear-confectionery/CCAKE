#ifndef KERNEL_H
#define KERNEL_H

#include <cmath>

#include <Cabana_Core.hpp>

namespace ccake{

/// @brief SPH kernel family, selected at runtime via the `kernel_type` setting.
/// Integer values are the `settings::kernel_id` codes set in input.cpp.
enum class KernelType : int {
  CUBIC       = 0,  ///< cubic B-spline (M4), support 2h  [default, legacy behaviour]
  QUINTIC     = 1,  ///< quintic B-spline (M6), support 3h
  WENDLAND_C2 = 2,  ///< Wendland C2, support 2h
  WENDLAND_C4 = 3,  ///< Wendland C4, support 2h
  WENDLAND_C6 = 4   ///< Wendland C6, support 2h
};

/// @class SPHkernel
/// @brief A class containing static functions for the SPH kernels.
/// @details This class provides methods for calculating the SPH kernel and
/// its gradient, as well as calculating the distance between two points.
/// The kernel family is chosen at call time by the @p kt argument (defaults to
/// the legacy cubic spline so untouched call sites are bit-identical).
/// @tparam D The dimension of the space.
template<unsigned int D>
class SPHkernel
{
private:
  static const double knorm, knorm2, kgrad, kgrad2;
  /// dimensionless shape f(q) and its derivative for the NON-cubic kernels
  /// (W = norm_C/h^D * f(q), q=r/h; D-aware for the Wendland 1D vs 2D/3D forms).
  KOKKOS_INLINE_FUNCTION static double shape_f(double q, KernelType kt);
  KOKKOS_INLINE_FUNCTION static double shape_dfdq(double q, KernelType kt);
  /// normalization C_D so that \int W d^Dx = 1 for the non-cubic kernels.
  KOKKOS_INLINE_FUNCTION static double norm_C(KernelType kt);

public:

  /// @brief Kernel support radius in units of h (the cutoff is support()*h):
  /// 2 for cubic and all Wendland kernels, 3 for the quintic spline.
  KOKKOS_INLINE_FUNCTION static double support(KernelType kt = KernelType::CUBIC)
  {
    return (kt == KernelType::QUINTIC) ? 3.0 : 2.0;
  }

  KOKKOS_FUNCTION static double kernel(double distance, double hT,
                                       KernelType kt = KernelType::CUBIC);
  KOKKOS_FUNCTION static double kernel(const double* r1, const double* r2,
                                       double hT, double hEta,
                                       KernelType kt = KernelType::CUBIC);
  /// @brief Partial derivative of the kernel w.r.t. the smoothing length, dW/dh.
  /// Needed for the grad-h (Springel-Hernquist/Price) Omega correction in adaptive-h.
  KOKKOS_FUNCTION static double dWdh(double distance, double h,
                                     KernelType kt = KernelType::CUBIC);
  KOKKOS_FUNCTION static void gradKernel(double const* rel_dist, double r,
                                         double h, double* grad,
                                         KernelType kt = KernelType::CUBIC);
  KOKKOS_FUNCTION static void gradKernel(const double* rel_sep,
                                         const double* pos_a, const double* pos_b,
                                         double hT, double hEta, double* grad,
                                         KernelType kt = KernelType::CUBIC);
  KOKKOS_FUNCTION static double distance(const double* r1, const double* r2);
};

}
#endif
