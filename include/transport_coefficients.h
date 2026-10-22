#ifndef TRANSPORT_COEFFICIENTS_H
#define TRANSPORT_COEFFICIENTS_H

#include <cmath>
#include <string>

#include <Cabana_Core.hpp>

#include "eos.h"
#include "settings.h"
#include "thermodynamic_info.h"
#include "constants.h"
#include "matrix.h"
#include "utilities.h"

namespace ccake{
namespace transport_coefficients{
  enum {
    SHEAR_DEFAULT,
    SHEAR_CONST,
    SHEAR_JAKI,
    SHEAR_LINEAR_MUS,
    SHEAR_PIECEWISE_MUB,   ///< BES-Bayesian piecewise-linear eta/s(muB), nodes at muB = 0, 0.2, 0.4 GeV
    SHEAR_INTERPOLATOR
  };
  enum {
    TAU_PI_SHEAR_DEFAULT,
    TAU_PI_SHEAR_MINVAL,
    TAU_PI_SHEAR_GUBSER,
  };
  enum {
    ZETA_DEFAULT,
    ZETA_CONSTANT,
    ZETA_DNMR,
    ZETA_INTERPOLATE,
    ZETA_CS2_DEPENDENT
  };
  enum {
    TAU_PI_BULK_DEFAULT,
    TAU_PI_BULK_DNMR
  };
  enum {
    KAPPA_DEFAULT,
    KAPPA_DNMR,
    KAPPA_B_COTH,   ///< baryon-only coth-form kappa_B (Eq. 16), needs B-only charges
    KAPPA_CONST_NB_OVER_T, ///< kappa_ij = const_ij * sqrt(|rho_i rho_j|)/T (density-weighted constant)
    KAPPA_B_HOLO_TABLE ///< baryon-only kappa_B = T*sigma_B(T,muB) from a tabulated holographic sigma_B
  };
  enum {
    TAU_Q_DEFAULT,      ///< tau_q = 0.2/T (constant_over_T)
    TAU_Q_CB_OVER_T,    ///< tau_q = C_B/T (couples relaxation to C_B; kappaB_coth legacy)
    TAU_Q_DNMR,
    TAU_Q_CATTANEO, ///< tau_q = D_B/v_0^2, D_B = kappa_BB/(T chi_BB) [Kapusta & Young, PRC 90, 044902]
    TAU_Q_KAPPA_OVER_N ///< tau_q = b_q * kappa_BB / n_B
  };
  enum {
    DELTA_B_PIPI_DEFAULT,
    DELTA_B_PIPI_IS,
    DELTA_B_PIPI_DISABLED
  };
  enum {
    LAMBDA_B_PIPI_DEFAULT,
    LAMBDA_B_PIPI_DISABLED
  };
  enum {
    PHI1_DEFAULT,
    PHI1_DISABLED
  };
  enum {
    PHI3_DEFAULT,
    PHI3_DISABLED
  };
  enum {
    DELTA_S_PIPI_DEFAULT,
    DELTA_S_PIPI_IS,
    DELTA_S_PIPI_DISABLED
  };
  enum {
    TAU_PIPI_DEFAULT,
    TAU_PIPI_DISABLED
  };
  enum {
    LAMBDA_S_PIPI_DEFAULT,
    LAMBDA_S_PIPI_DISABLED
  };
  enum {
    PHI6_DEFAULT,
    PHI6_DISABLED
  };
  enum {
    PHI7_DEFAULT,
    PHI7_DISABLED
  };
    

  struct parameters
  {
    int shear_mode, zeta_mode, shear_relaxation_mode, bulk_relaxation_mode, diffusion_mode, diffusion_relaxation_mode;
    int delta_pipi_mode, tau_pipi_mode, lambda_piPi_mode, phi6_mode, phi7_mode;
    int delta_PiPi_mode, lambda_Pipi_mode, phi1_mode, phi3_mode;
    double constant_eta_over_s, constant_zeta_over_s;
    double shear_eta_0, shear_eta_2, shear_eta_4; ///< eta/s at muB = 0, 0.2, 0.4 GeV (piecewise_muB shear mode)
    double cs2_dependent_zeta_A, cs2_dependent_zeta_p, bulk_tau_pi_coeff;
    double C_B;   ///< coefficient for coth-form baryon diffusion kappa_B (Eq. 16)
    ///@{ holographic sigma_B(T, muB) table (diffusion mode holographic_table, kappa_B = scale T sigma_B),
    ///   indexed [T][muB], natural units; num_temperatures == 0 means no table.
    Kokkos::View<double**, MemorySpace> holo_sigmaB;
    int    holo_sigmaB_num_temperatures = 0;
    int    holo_sigmaB_num_potentials   = 0;
    double holo_sigmaB_temperature_min = 0.0, holo_sigmaB_temperature_step = 0.0;
    double holo_sigmaB_potential_min   = 0.0, holo_sigmaB_potential_step   = 0.0;
    double holo_diffusion_scale = 1.0; ///< dimensionless overall multiplier on kappa_B
    ///@}
    ///@{ cattaneo mode: tau_q = D_B/diffusion_wave_speed_sq (<= 0: local cs2); bounds in units of 1/T.
    double diffusion_wave_speed_sq = 0.0;
    double tau_q_min_over_T = 0.1;
    double tau_q_max_over_T = 10.0;
    ///< kappa_over_n relaxation mode: tau_q = b_q * kappa_BB/n_B.
    double b_q = 20.0;
    ///< enforce the Kapusta-Young causality floor tau_q >= D_B (default true)
    bool tau_q_causality_floor = true;
    ///@}
    std::array<std::array<double, 3>, 3> constant_kappa_over_T2,
                                         tauq_matrix;
    bool modulate_zeta_with_tanh;
    bool critical_scaling_bulk;
    /// use the table's local chiBB_nc baseline as xi_0 (true) or the global constant (false)
    bool critical_scaling_local_xi0 = false;
    /// hard cap on the bulk critical enhancement; <= 0 disables the cap
    double critical_scaling_max = 10.0;
    bool critical_scaling_diffusion;
    double critical_point_T;            ///< Critical point T in natural units
    double critical_point_muB;          ///< Critical point muB in natural units
    double critical_gaussian_width_T;   ///< Gaussian width in T in natural units
    double critical_gaussian_width_muB; ///< Gaussian width in muB in natural units
    /// 0 = legacy Gaussian suppression of kappa_B, 1 = xi_power (kappa_B and tau_q enhanced)
    int critical_scaling_diffusion_mode = 0;
    double critical_scaling_kappa_power  = 1.0;  ///< kappa_B ~ xi^p in xi_power mode
    double critical_scaling_tauq_power   = 2.0;  ///< tau_q ~ xi^p in xi_power mode
  };

  KOKKOS_INLINE_FUNCTION
  parameters setup_parameters(std::shared_ptr<Settings> settingsPtr);

  /// Host-only (file I/O): load the holographic sigma_B table into \p params after
  /// setup_parameters(). No-op unless the diffusion mode is "holographic_table".
  inline void load_holographic_sigmaB_table(parameters & params,
                                            std::shared_ptr<Settings> settingsPtr);
  KOKKOS_INLINE_FUNCTION
  double eta(const double* thermo, parameters params);
  KOKKOS_INLINE_FUNCTION
  double tau_pi(const double* thermo, parameters params);
  KOKKOS_INLINE_FUNCTION
  double zeta(const double* thermo, parameters params);
  KOKKOS_INLINE_FUNCTION
  double tau_Pi(const double* thermo, parameters params);
  KOKKOS_INLINE_FUNCTION
  double delta_PiPi(const double *therm, parameters params);
  KOKKOS_INLINE_FUNCTION
  double lambda_Pipi(const double *therm, parameters params);
  KOKKOS_INLINE_FUNCTION
  double phi1(const double *therm, parameters params);
  KOKKOS_INLINE_FUNCTION
  double phi3(const double *therm, parameters params);
  KOKKOS_INLINE_FUNCTION
  double delta_pipi(const double *therm, parameters params);
  KOKKOS_INLINE_FUNCTION
  double tau_pipi(const double *therm, parameters params);
  KOKKOS_INLINE_FUNCTION
  double lambda_piPi(const double *therm, parameters params);
  KOKKOS_INLINE_FUNCTION
  double phi6(const double *therm, parameters params);
  KOKKOS_INLINE_FUNCTION
  double phi7(const double *therm, parameters params);
  //critical functions
  KOKKOS_INLINE_FUNCTION
  double critical_scaling_factor_bulk(const double *therm, parameters params);
  KOKKOS_INLINE_FUNCTION
  double gaussian_suppression_diffusion(const double *therm, parameters params);
  KOKKOS_INLINE_FUNCTION
  double critical_xi_power_factor(const double *therm, parameters params, double power);
  KOKKOS_INLINE_FUNCTION
  double critical_factor_kappa(const double *therm, parameters params);
  KOKKOS_INLINE_FUNCTION
  double critical_factor_tauq(const double *therm, parameters params);


  // Chris' defaults
  KOKKOS_INLINE_FUNCTION
  double default_eta(const double *therm);
  KOKKOS_INLINE_FUNCTION
  double default_tau_pi(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  double default_zeta(const double *therm);
  KOKKOS_INLINE_FUNCTION
  double default_tau_Pi(const double *therm, const parameters params);

  //DNMR terms
  KOKKOS_INLINE_FUNCTION
  double delta_PiPi_DEFAULT(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  double delta_PiPi_IS(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  double lambda_Pipi_DEFAULT(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  double phi1_DEFAULT(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION  
  double phi3_DEFAULT(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  double delta_pipi_DEFAULT(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  double delta_pipi_IS(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  double tau_pipi_DEFAULT(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  double lambda_piPi_DEFAULT(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  double phi6_DEFAULT(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION    
  double phi7_DEFAULT(const double *therm, const parameters params);

  // diffusion
  KOKKOS_INLINE_FUNCTION
  Matrix<double, 3, 3> kappa(const double* thermo, parameters params);
  KOKKOS_INLINE_FUNCTION
  Matrix<double, 3, 3> tauq(const double* thermo, parameters params);
  KOKKOS_INLINE_FUNCTION
  Matrix<double, 3, 3> default_kappa(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  Matrix<double, 3, 3> default_tauq(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  Matrix<double, 3, 3> kappaB_coth(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  Matrix<double, 3, 3> const_nB_over_T_kappa(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  Matrix<double, 3, 3> kappaB_holo_table(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  /// @brief Raises tau_q to the causality bound tau_q >= D_B = kappa_BB/(T chi_BB), if
  /// tau_q_causality_floor (cattaneo and kappa_over_n enforce it themselves).
  double tau_q_with_causality_floor(double tau_q, const double *therm,
                                    const parameters params);

  Matrix<double, 3, 3> tauq_cattaneo(const double *therm, const parameters params);
  Matrix<double, 3, 3> tauq_kappa_over_n(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  Matrix<double, 3, 3> tauqB_only(const double *therm, const parameters params);

  // other parmaterizations (not organized yet)
  KOKKOS_INLINE_FUNCTION
  double constEta(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  double JakiParam(const double *therm);
  KOKKOS_INLINE_FUNCTION
  double LinearMusParam(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  double PiecewiseMuBParam(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  double InterpolantWrapper(const double *therm);

  KOKKOS_INLINE_FUNCTION
  double tau_piGubser(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  double tau_piMinval(const double *therm, const parameters params);

  KOKKOS_INLINE_FUNCTION
  double constZeta(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  double zeta_DNMR_LeadingMass(const double *therm);
  KOKKOS_INLINE_FUNCTION
  double cs2_dependent_zeta(const double *therm, const parameters params);

  KOKKOS_INLINE_FUNCTION
  double tau_Pi_DNMR15(const double *therm, const parameters params);

}}


#endif
