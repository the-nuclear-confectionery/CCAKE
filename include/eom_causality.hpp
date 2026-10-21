#ifndef EOM_CAUSALITY_HPP
#define EOM_CAUSALITY_HPP

/// @file eom_causality.hpp
/// @brief Nonlinear causality check (Bemfica et al., arXiv:2005.11632), shared by the Milne and
/// Cartesian EoM: necessary plus sufficient conditions, hence causal / undetermined /
/// acausal. CCAKE is mostly-minus and stores pi^{mu nu} contravariant; the paper is
/// mostly-plus, so Lambda_paper = -Lambda_CCAKE (handled in shear_eigenvalues).

#include <memory>

#include <Kokkos_Core.hpp>

#include "particle.h"
#include "system_state.h"

namespace ccake {
namespace causality {

/// Per-particle verdict codes (stored in hydro_info::causality).
enum verdict : int
{
  ERROR        = -100, ///< inputs degenerate; verdict not computable
  ACAUSAL      = -1,
  UNDETERMINED = 0,    ///< necessary conditions hold, sufficient ones do not
  CAUSAL       = 1
};

namespace detail
{
  constexpr double TINY = 1.0e-30;
  constexpr double TWO_PI_OVER_3 = 2.0943951023931954;
}

/// @brief Inputs of the verdict, in the paper's conventions.
struct Inputs
{
  double Lam[3];   ///< eigenvalues of pi^mu_nu, ascending, sum = 0
  double cs2;      ///< speed of sound squared
  double e, p, Pi;
  double eta, zeta, tau_pi, tau_Pi;
  double delta_PiPi, lambda_Pipi;   ///< bulk equation:  delta_PiPi, lambda_Ppi
  double delta_pipi, tau_pipi, lambda_piPi;  ///< shear equation
};

//==============================================================================
/// @brief Eigenvalues Lambda_1 <= Lambda_2 <= Lambda_3 of the mixed shear tensor.
/// pi^mu_nu is traceless and annihilates u^mu, so they are the roots of the cubic
/// lambda^3 - p lambda - q = 0 (p = tr(pi^2)/2, q = tr(pi^3)/3), solved in closed form.
/// q flips sign between the metric conventions.
/// @param shv contravariant pi^{mu nu} (hydro.shv); @param t2 tau^2; @param Lam output, ascending.
template<class MatrixT>
KOKKOS_INLINE_FUNCTION
void shear_eigenvalues(const MatrixT& shv, double t2, double Lam[3])
{
  MatrixT mixed = shv;
  mixed.make_covariant(1, t2);          // pi^{mu nu} -> pi^mu_nu (lower column index)

  double M[4][4];
  for (int idir = 0; idir < 4; ++idir)
    for (int jdir = 0; jdir < 4; ++jdir)
      M[idir][jdir] = mixed(idir, jdir);

  double M2[4][4];
  for (int idir = 0; idir < 4; ++idir)
    for (int jdir = 0; jdir < 4; ++jdir)
    {
      double sum = 0.0;
      for (int k = 0; k < 4; ++k) sum += M[idir][k]*M[k][jdir];
      M2[idir][jdir] = sum;
    }

  double trace2 = 0.0, trace3 = 0.0;
  for (int idir = 0; idir < 4; ++idir) trace2 += M2[idir][idir];
  for (int idir = 0; idir < 4; ++idir)
    for (int k = 0; k < 4; ++k) trace3 += M2[idir][k]*M[k][idir];

  const double p =  0.5*trace2;
  const double q = -trace3/3.0;         // sign flip: CCAKE mostly-minus -> paper mostly-plus

  if ( p <= detail::TINY )              // vanishing shear stress
  {
    Lam[0] = Lam[1] = Lam[2] = 0.0;
    return;
  }

  double arg = (3.0*q)/(2.0*p)*Kokkos::sqrt(3.0/p);
  arg = Kokkos::fmax(-1.0, Kokkos::fmin(1.0, arg));   // guard round-off at |arg| = 1
  const double theta  = Kokkos::acos(arg)/3.0;
  const double radius = 2.0*Kokkos::sqrt(p/3.0);

  double root[3];
  for (int k = 0; k < 3; ++k)
    root[k] = radius*Kokkos::cos(theta - detail::TWO_PI_OVER_3*k);

  // sort ascending (3 elements; a sorting network is cheaper than a loop)
  if (root[0] > root[1]) { const double s = root[0]; root[0] = root[1]; root[1] = s; }
  if (root[1] > root[2]) { const double s = root[1]; root[1] = root[2]; root[2] = s; }
  if (root[0] > root[1]) { const double s = root[0]; root[0] = root[1]; root[1] = s; }

  Lam[0] = root[0]; Lam[1] = root[1]; Lam[2] = root[2];
}

//==============================================================================
/// @brief Necessary (4a-4f) and sufficient (5a-5h) conditions of arXiv:2005.11632, at zero density
/// (cs2 stands in for P_eps). @return CAUSAL / UNDETERMINED / ACAUSAL.
KOKKOS_INLINE_FUNCTION
int classify(const Inputs& in)
{
  const double Lambda_1 = in.Lam[0], Lambda_2 = in.Lam[1], Lambda_3 = in.Lam[2];
  const double e = in.e, p = in.p, Pi = in.Pi, cs2 = in.cs2;
  const double eta = in.eta, zeta = in.zeta;
  const double tau_pi = in.tau_pi, tau_Pi = in.tau_Pi;
  const double tau_pipi = in.tau_pipi, delta_PiPi = in.delta_PiPi;
  const double lambda_Pipi = in.lambda_Pipi, delta_pipi = in.delta_pipi;
  const double lambda_piPi = in.lambda_piPi;

  if ( tau_pi <= detail::TINY || tau_Pi <= detail::TINY ) return ERROR;

  const double shear_drive = 2.0*eta + lambda_piPi*Pi;   // 2*eta_eff

  bool necessary
      = (tau_Pi>=0) && (tau_pi>=0) && (eta>=0) && (zeta>=0)
        && (tau_pipi>=0) && (delta_PiPi>=0) && (lambda_Pipi>=0)
        && (delta_pipi>=0) && (lambda_piPi>=0) && (cs2>=0)
        && (e>=0) && (p>=0) && (e+p+Pi>0);

  const bool condition4a = shear_drive - tau_pi*Kokkos::fabs(Lambda_1)/2.0 >= 0;
  const bool condition4b = e + p + Pi - shear_drive/(2.0*tau_pi)
                           - (tau_pipi/(4.0*tau_pi))*Lambda_3 >= 0;
  necessary = necessary && condition4a && condition4b;

  for ( int d = 0; d < 3; ++d )
  {
    const double Lambda_d = in.Lam[d];
    const double common = shear_drive/(2.0*tau_pi)
                        + (tau_pipi/(2.0*tau_pi))*Lambda_d
                        + (shear_drive + (6.0*delta_pipi - tau_pipi)*Lambda_d)/(6.0*tau_pi)
                        + (zeta + delta_PiPi*Pi + lambda_Pipi*Lambda_d)/tau_Pi
                        + (e + p + Pi + Lambda_d)*cs2;
    const bool condition4e = common >= 0;
    const bool condition4f = (e + p + Pi + Lambda_d) - common >= 0;
    necessary = necessary && condition4e && condition4f;
  }

  for ( int a = 0; a < 3; ++a )
    for ( int d = 0; d < 3; ++d )
    {
      if ( a == d ) continue;
      const double Lambda_a = in.Lam[a], Lambda_d = in.Lam[d];
      const bool condition4c = shear_drive/(2.0*tau_pi)
                             + (tau_pipi/(4.0*tau_pi))*(Lambda_a + Lambda_d) >= 0;
      const bool condition4d = e + p + Pi + Lambda_a - shear_drive/(2.0*tau_pi)
                             - (tau_pipi/(4.0*tau_pi))*(Lambda_d + Lambda_a) >= 0;
      necessary = necessary && condition4c && condition4d;
    }

  if ( !necessary ) return ACAUSAL;

  const double abs_L1 = Kokkos::fabs(Lambda_1);
  const double rest   = e + p + Pi - abs_L1;
  const double margin = rest - shear_drive/(2.0*tau_pi)
                        - (tau_pipi/(2.0*tau_pi))*Lambda_3;

  const bool condition5b = shear_drive - tau_pipi*abs_L1 > 0;
  const bool condition5c = tau_pipi <= 6.0*delta_pipi;
  const bool condition5d = lambda_Pipi/tau_Pi + cs2 - tau_pipi/(12.0*tau_pi) >= 0;

  // 5a (margin >= 0) is required strictly, since 5e divides by margin and 5h by rest
  if ( !(condition5b && condition5c && condition5d)
       || margin <= detail::TINY || rest <= detail::TINY ) return UNDETERMINED;

  const double gap  = Lambda_3 + abs_L1;
  const double curv = ((12.0*delta_pipi - tau_pipi)/(12.0*tau_pi))
                      * (lambda_Pipi/tau_Pi + cs2 - tau_pipi/(12.0*tau_pi))*gap*gap;

  const bool condition5e =
      (4.0*eta + 2.0*lambda_piPi*Pi + (3.0*delta_pipi + tau_pipi)*Lambda_3)/(3.0*tau_pi)
      + (zeta + delta_PiPi*Pi + lambda_Pipi*Lambda_3)/tau_Pi
      + abs_L1 + Lambda_3*cs2 + curv/margin
      <= (e + p + Pi)*(1.0 - cs2);

  const bool condition5f =
      (shear_drive + (tau_pipi - 6.0*delta_pipi)*abs_L1)/(6.0*tau_pi)
      + (zeta + delta_PiPi*Pi - lambda_Pipi*abs_L1)/tau_Pi
      + rest*cs2 >= 0;

  // condition5b above is exactly denom5g > 0, so this division is safe
  const double denom5g = shear_drive/(2.0*tau_pi) - (tau_pipi/(2.0*tau_pi))*abs_L1;
  const bool condition5g = 1.0 >= curv/(denom5g*denom5g);

  const bool condition5h =
      (4.0*eta + 2.0*lambda_piPi*Pi - (3.0*delta_pipi + tau_pipi)*abs_L1)/(3.0*tau_pi)
      + (zeta + delta_PiPi*Pi - lambda_Pipi*abs_L1)/tau_Pi
      + rest*cs2
      >= ((e + p + Pi + Lambda_2)*(e + p + Pi + Lambda_3)) / (3.0*rest)
         * (1.0 + 2.0*(shear_drive/(2.0*tau_pi) + (tau_pipi/(2.0*tau_pi))*Lambda_3)
                  / rest);

  return (condition5e && condition5f && condition5g && condition5h) ? CAUSAL : UNDETERMINED;
}

//==============================================================================
/// @brief Classify every SPH particle and store the verdict in hydro_info::causality.
/// MatrixT is milne::Matrix<double,4,4> or cartesian::Matrix<double,4,4>.
template <class MatrixT, unsigned int D>
void classify_particles(std::shared_ptr<SystemState<D>> sysPtr)
{
  const double t  = sysPtr->t;
  const double t2 = t*t;

  CREATE_VIEW(device_, sysPtr->cabana_particles);
  auto host_simd_policy
      = Cabana::SimdPolicy<VECTOR_LENGTH, Kokkos::OpenMP>(0, sysPtr->cabana_particles.size());

  Kokkos::fence(); // ensure device kernels have completed before the host-side check

  auto causality_check = [=](const int is, const int ia)
  {
    Inputs in;

    in.e   = device_thermo.access(is, ia, thermo_info::e);
    in.p   = device_thermo.access(is, ia, thermo_info::p);
    in.cs2 = device_thermo.access(is, ia, thermo_info::cs2);
    in.Pi  = device_hydro_scalar.access(is, ia, hydro_info::bulk);

    in.eta         = device_hydro_scalar.access(is, ia, hydro_info::eta_pi);
    in.zeta        = device_hydro_scalar.access(is, ia, hydro_info::zeta_Pi);
    in.tau_pi      = device_hydro_scalar.access(is, ia, hydro_info::tau_pi);
    in.tau_Pi      = device_hydro_scalar.access(is, ia, hydro_info::tau_Pi);
    in.tau_pipi    = device_hydro_scalar.access(is, ia, hydro_info::tau_pipi);
    in.delta_PiPi  = device_hydro_scalar.access(is, ia, hydro_info::delta_PiPi);
    in.lambda_Pipi = device_hydro_scalar.access(is, ia, hydro_info::lambda_Pipi);
    in.delta_pipi  = device_hydro_scalar.access(is, ia, hydro_info::delta_pipi);
    in.lambda_piPi = device_hydro_scalar.access(is, ia, hydro_info::lambda_piPi);

    MatrixT shv;
    for ( int idir = 0; idir < 4; ++idir )
      for ( int jdir = 0; jdir < 4; ++jdir )
        shv(idir, jdir)
            = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir);

    shear_eigenvalues(shv, t2, in.Lam);

    device_hydro_scalar.access(is, ia, hydro_info::causality) = classify(in);
  };

  Cabana::simd_parallel_for(host_simd_policy, causality_check, "check_causality_conditions");
  Kokkos::fence();
}

} // namespace causality
} // namespace ccake

#endif
