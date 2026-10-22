#ifndef TRANSPORT_COEFFICIENTS_CPP
#define TRANSPORT_COEFFICIENTS_CPP
#include "transport_coefficients.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "formatted_output.h"

namespace ccake{
namespace transport_coefficients{

KOKKOS_INLINE_FUNCTION
parameters setup_parameters(std::shared_ptr<Settings> settingsPtr)
{
  std::string etaMode = settingsPtr->etaMode;
  parameters params;
  if (etaMode == "default")
  {
    params.shear_mode = SHEAR_DEFAULT;
  } else if (etaMode == "JakiParam")
  {
    params.shear_mode = SHEAR_JAKI;
  }
  else if ( etaMode == "constant" )
  {
      params.shear_mode = SHEAR_CONST;
  }
  else if ( etaMode == "LinearMus" )
  {
    params.shear_mode = SHEAR_LINEAR_MUS;
  }
  else if ( etaMode == "piecewise_muB" )
  {
    params.shear_mode = SHEAR_PIECEWISE_MUB;
  }
  else if ( etaMode == "interpolate" )
  {
    params.shear_mode = SHEAR_INTERPOLATOR;
  }

  std::string shearRelaxMode = settingsPtr->shearRelaxMode;
  if (shearRelaxMode == "default")
  {
    params.shear_relaxation_mode = TAU_PI_SHEAR_DEFAULT;
  }
  else if (shearRelaxMode == "minVal")
  {
    params.shear_relaxation_mode = TAU_PI_SHEAR_MINVAL;
  }
  else if (shearRelaxMode == "Gubser")
  {
    params.shear_relaxation_mode = TAU_PI_SHEAR_GUBSER;
  }

  std::string zetaMode = settingsPtr->zetaMode;
  if (zetaMode == "default")
  {
    params.zeta_mode = ZETA_DEFAULT;
  }
  else if (zetaMode == "constant")
  {
    params.zeta_mode = ZETA_CONSTANT;
  }
  else if (zetaMode == "DNMR")
  {
    params.zeta_mode = ZETA_DNMR;
  }
  else if (zetaMode == "interpolate"){
    params.zeta_mode = ZETA_INTERPOLATE;
  }
  else if (zetaMode == "cs2_dependent"){
    params.zeta_mode = ZETA_CS2_DEPENDENT;
  }
  params.critical_scaling_local_xi0 = settingsPtr->critical_scaling_local_xi0;
  params.critical_scaling_max       = settingsPtr->critical_scaling_max;
  if (settingsPtr->critical_scaling_bulk) {
    params.critical_scaling_bulk = true;
  } else {
    params.critical_scaling_bulk = false;
  }

  std::string bulkRelaxMode = settingsPtr->bulkRelaxMode;
  if (bulkRelaxMode == "default")
  {
    params.bulk_relaxation_mode = TAU_PI_BULK_DEFAULT;
  }
  else if(bulkRelaxMode == "DNMR")
  {
    params.bulk_relaxation_mode = TAU_PI_BULK_DNMR;
  }

  std::string diffusionMode = settingsPtr->diffusionMode;
  if (diffusionMode == "constant_over_T2")
  {
    params.diffusion_mode = KAPPA_DEFAULT;
  }
  else if (diffusionMode == "DNMR")
  {
    params.diffusion_mode = KAPPA_DNMR;
  }
  else if (diffusionMode == "kappaB_coth")
  {
    params.diffusion_mode = KAPPA_B_COTH;
  }
  else if (diffusionMode == "constant_nB_over_T")
  {
    params.diffusion_mode = KAPPA_CONST_NB_OVER_T;
  }
  else if (diffusionMode == "holographic_table")
  {
    params.diffusion_mode = KAPPA_B_HOLO_TABLE;
  }
  else
  {
    std::cout << "WARNING: Unknown diffusion mode: " << diffusionMode << std::endl;
    std::cout << "Defaulting to constant_over_T2." << std::endl;
    params.diffusion_mode = KAPPA_DEFAULT; // default to KAPPA_DEFAULT if not specified
  }

  std::string diffusionRelaxMode = settingsPtr->diffusionRelaxMode;
  if (diffusionRelaxMode == "constant_over_T")
  {
    params.diffusion_relaxation_mode = TAU_Q_DEFAULT;      // tau_q = 0.2/T
  }
  else if (diffusionRelaxMode == "cb_over_T")
  {
    params.diffusion_relaxation_mode = TAU_Q_CB_OVER_T;    // tau_q = C_B/T
  }
  else if (diffusionRelaxMode == "kappa_over_n")
  {
    params.diffusion_relaxation_mode = TAU_Q_KAPPA_OVER_N;
  }
  else if (diffusionRelaxMode == "cattaneo" || diffusionRelaxMode == "nernst_einstein")
  {
    params.diffusion_relaxation_mode = TAU_Q_CATTANEO;
  }
  else{
    std::cout << "WARNING: Unknown diffusion relaxation mode: " << diffusionRelaxMode << std::endl
              << "Defaulting to constant_over_T (tau_q = 0.2/T)." << std::endl;
    params.diffusion_relaxation_mode = TAU_Q_DEFAULT; // default to TAU_Q_DEFAULT if not specified

  }
  std::string delta_pipi_mode = settingsPtr->delta_pipi_mode;
  if (delta_pipi_mode == "default")
  {
    params.delta_pipi_mode = DELTA_S_PIPI_DEFAULT;
  }
  else if (delta_pipi_mode == "israel-stewart")
  {
    params.delta_pipi_mode = DELTA_S_PIPI_IS;
  }
  else if (delta_pipi_mode == "disabled")
  {
    params.delta_pipi_mode = DELTA_S_PIPI_DISABLED;
  }
  else
  {
    std::cout << "WARNING: Unknown delta_pipi mode: " << delta_pipi_mode << std::endl;
    std::cout << "Defaulting to delta_pipi_DEFAULT." << std::endl;
    params.delta_pipi_mode = DELTA_S_PIPI_DEFAULT; // default to delta_pipi_DEFAULT if not specified
  }
  std::string tau_pipi_mode = settingsPtr->tau_pipi_mode;
  if (tau_pipi_mode == "default")
  {
    params.tau_pipi_mode = TAU_PIPI_DEFAULT;
  }
  else if (tau_pipi_mode == "disabled")
  {
    params.tau_pipi_mode = TAU_PIPI_DISABLED;
  }
  else
  {
    std::cout << "WARNING: Unknown tau_pipi mode: " << tau_pipi_mode << std::endl;
    std::cout << "Defaulting to tau_pipi_DEFAULT." << std::endl;
    params.tau_pipi_mode = TAU_PIPI_DEFAULT; // default to tau_pipi_DEFAULT if not specified
  }
  std::string lambda_piPi_mode = settingsPtr->lambda_piPi_mode;
  if (lambda_piPi_mode == "default")
  {
    params.lambda_piPi_mode = LAMBDA_S_PIPI_DEFAULT;
  
  }
  else if (lambda_piPi_mode == "disabled")
  {
    params.lambda_piPi_mode = LAMBDA_S_PIPI_DISABLED;
  }
  else
  {
    std::cout << "WARNING: Unknown lambda_piPi mode: " << lambda_piPi_mode << std::endl;
    std::cout << "Defaulting to lambda_piPi_DEFAULT." << std::endl;
    params.lambda_piPi_mode = LAMBDA_S_PIPI_DEFAULT; // default to lambda_PiPi_DEFAULT if not specified
  }
  std::string phi6_mode = settingsPtr->phi6_mode;
  if (phi6_mode == "default")
  {
    params.phi6_mode = PHI6_DEFAULT;
  }
  else if (phi6_mode == "disabled")
  {
    params.phi6_mode = PHI6_DISABLED;
  }
  else
  {
    std::cout << "WARNING: Unknown phi6 mode: " << phi6_mode << std::endl;
    std::cout << "Defaulting to phi6_DEFAULT." << std::endl;
    params.phi6_mode = PHI6_DEFAULT; // default to phi6_DEFAULT if not specified
  }
  std::string phi7_mode = settingsPtr->phi7_mode;
  if (phi7_mode == "default")
  {
    params.phi7_mode = PHI7_DEFAULT;
  }
  else if (phi7_mode == "disabled")
  {
    params.phi7_mode = PHI7_DISABLED;
  }
  else
  {
    std::cout << "WARNING: Unknown phi7 mode: " << phi7_mode << std::endl;
    std::cout << "Defaulting to phi7_DEFAULT." << std::endl;
    params.phi7_mode = PHI7_DEFAULT; // default to phi7_DEFAULT if not specified
  }

  std::string delta_PiPi_mode = settingsPtr->delta_PiPi_mode;
  if (delta_PiPi_mode == "default")
  {
    params.delta_PiPi_mode = DELTA_B_PIPI_DEFAULT;
  }
  else if (delta_PiPi_mode == "israel-stewart")
  {
    params.delta_PiPi_mode = DELTA_B_PIPI_IS;
  }
  else if (delta_PiPi_mode == "disabled")
  {
    params.delta_PiPi_mode = DELTA_B_PIPI_DISABLED;
  }
  else
  {
    std::cout << "WARNING: Unknown delta_PiPi mode: " << delta_PiPi_mode << std::endl;
    std::cout << "Defaulting to delta_PiPi_DEFAULT." << std::endl;
    params.delta_PiPi_mode = DELTA_B_PIPI_DEFAULT; // default to delta_PiPi_DEFAULT if not specified
  }
  std::string lambda_Pipi_mode = settingsPtr->lambda_Pipi_mode;
  if (lambda_Pipi_mode == "default")
  {
    params.lambda_Pipi_mode = LAMBDA_B_PIPI_DEFAULT;
  }
  else if (lambda_Pipi_mode == "disabled")
  {
    params.lambda_Pipi_mode = LAMBDA_B_PIPI_DISABLED;
  }
  else
  {
    std::cout << "WARNING: Unknown lambda_Pipi mode: " << lambda_Pipi_mode << std::endl;
    std::cout << "Defaulting to lambda_Pipi_DEFAULT." << std::endl;
    params.lambda_Pipi_mode = LAMBDA_B_PIPI_DEFAULT; // default to lambda_PiPi_DEFAULT if not specified
  }
  std::string phi1_mode = settingsPtr->phi1_mode;
  if (phi1_mode == "default")
  {
    params.phi1_mode = PHI1_DEFAULT;
  }
  else if (phi1_mode == "disabled")
  {
    params.phi1_mode = PHI1_DISABLED;
  }
  else
  {
    std::cout << "WARNING: Unknown phi1 mode: " << phi1_mode << std::endl;
    std::cout << "Defaulting to phi1_DEFAULT." << std::endl;
    params.phi1_mode = PHI1_DEFAULT; // default to phi1_DEFAULT if not specified
  }
  std::string phi3_mode = settingsPtr->phi3_mode;
  if (phi3_mode == "default")
  {
    params.phi3_mode = PHI3_DEFAULT;
  }
  else if (phi3_mode == "disabled")
  {
    params.phi3_mode = PHI3_DISABLED;
  }
  else
  {
    std::cout << "WARNING: Unknown phi3 mode: " << phi3_mode << std::endl;
    std::cout << "Defaulting to phi3_DEFAULT." << std::endl;
    params.phi3_mode = PHI3_DEFAULT; // default to phi3_DEFAULT if not specified
  }




  params.constant_eta_over_s = settingsPtr->constant_eta_over_s;
  params.shear_eta_0 = settingsPtr->shear_eta_0;
  params.shear_eta_2 = settingsPtr->shear_eta_2;
  params.shear_eta_4 = settingsPtr->shear_eta_4;
  params.constant_zeta_over_s = settingsPtr->constant_zeta_over_s;
  params.cs2_dependent_zeta_A = settingsPtr->cs2_dependent_zeta_A;
  params.cs2_dependent_zeta_p = settingsPtr->cs2_dependent_zeta_p;
  params.bulk_tau_pi_coeff = settingsPtr->bulk_tau_pi_coeff;
  params.modulate_zeta_with_tanh = settingsPtr->modulate_zeta_with_tanh;
  params.C_B = settingsPtr->C_B;
  params.holo_diffusion_scale = settingsPtr->holographic_diffusion_scale;
  params.diffusion_wave_speed_sq = settingsPtr->diffusion_wave_speed_sq;
  params.b_q                     = settingsPtr->b_q;
  params.tau_q_causality_floor   = settingsPtr->tau_q_causality_floor;
  params.tau_q_min_over_T        = settingsPtr->tau_q_min_over_T;
  params.tau_q_max_over_T        = settingsPtr->tau_q_max_over_T;
  params.constant_kappa_over_T2 = settingsPtr->constant_kappa_over_T2;
  params.critical_scaling_diffusion  = settingsPtr->critical_scaling_diffusion;
  params.critical_point_T            = settingsPtr->critical_point_T   / constants::hbarc_MeVfm;
  params.critical_point_muB          = settingsPtr->critical_point_muB / constants::hbarc_MeVfm;
  params.critical_gaussian_width_T   = settingsPtr->critical_gaussian_width_T   / constants::hbarc_MeVfm;
  params.critical_gaussian_width_muB = settingsPtr->critical_gaussian_width_muB / constants::hbarc_MeVfm;
  params.critical_scaling_diffusion_mode = ( settingsPtr->critical_scaling_diffusion_mode == "xi_power" ) ? 1 : 0;
  params.critical_scaling_kappa_power    = settingsPtr->critical_scaling_kappa_power;
  params.critical_scaling_tauq_power     = settingsPtr->critical_scaling_tauq_power;

  return params;

}

KOKKOS_INLINE_FUNCTION
double eta(const double* thermo, parameters params)
{
  int eta_mode = params.shear_mode;
  switch (eta_mode)
  {
  case SHEAR_DEFAULT:
    return default_eta(thermo);
    break;
  case SHEAR_CONST:
    return constEta(thermo, params);
    break;
  case SHEAR_JAKI:
    return JakiParam(thermo);
    break;
  case SHEAR_LINEAR_MUS:
    return LinearMusParam(thermo, params);
    break;
  case SHEAR_PIECEWISE_MUB:
    return PiecewiseMuBParam(thermo, params);
    break;
  case SHEAR_INTERPOLATOR:
    return InterpolantWrapper(thermo);
    break;
  default:
    return default_eta(thermo);
    break;
  }
}

KOKKOS_INLINE_FUNCTION
double tau_pi(const double* thermo, parameters params)
{
  int tau_pi_mode = params.shear_relaxation_mode;
  switch (tau_pi_mode)
  {
  case TAU_PI_SHEAR_DEFAULT:
    return default_tau_pi(thermo, params);
    break;
  case TAU_PI_SHEAR_GUBSER:
    return tau_piGubser(thermo, params);
    break;
  case TAU_PI_SHEAR_MINVAL:
    return tau_piMinval(thermo, params);
    break;
  default:
    return default_tau_pi(thermo, params);
    break;
  }
}

KOKKOS_INLINE_FUNCTION
double zeta(const double* thermo, parameters params)
{
  int zeta_mode = params.zeta_mode;
  double zeta = 0.;
  switch (zeta_mode)
  {
  case ZETA_DEFAULT:
    zeta = default_zeta(thermo);
    break;
  case ZETA_CONSTANT:
    zeta = constZeta(thermo, params);
    break;
  case ZETA_DNMR:
    zeta = zeta_DNMR_LeadingMass(thermo);
    break;
  case ZETA_CS2_DEPENDENT:
    zeta = cs2_dependent_zeta(thermo, params);
    break;
  default:
    zeta = default_zeta(thermo);
    break;
  }
  //check for critical scaling
  if (params.critical_scaling_bulk) {
    //apply critical scaling factor
    zeta *= critical_scaling_factor_bulk(thermo, params);
  }
  return zeta;

}

KOKKOS_INLINE_FUNCTION
double tau_Pi(const double* thermo, parameters params)
{
  int tau_Pi_mode = params.bulk_relaxation_mode;
  switch (tau_Pi_mode)
  {
  case TAU_PI_BULK_DEFAULT:
    return default_tau_Pi(thermo,params);
    break;
  case TAU_PI_BULK_DNMR:
    return tau_Pi_DNMR15(thermo, params);   // tau_Pi = zeta / [15 (1/3-cs2)^2 w]
    break;
  default:
    return default_tau_Pi(thermo,params);
    break;
  }
}


KOKKOS_INLINE_FUNCTION
double delta_PiPi(const double* thermo, parameters params)
{
  int delta_PiPi_mode = params.delta_PiPi_mode;
  switch (delta_PiPi_mode)
  {
  case DELTA_B_PIPI_DEFAULT:
    return delta_PiPi_DEFAULT(thermo, params);
    break;
  case DELTA_B_PIPI_IS:
    return delta_PiPi_IS(thermo, params);
    break;
  case DELTA_B_PIPI_DISABLED:
    return 0.0; // No delta PiPi contribution
    break;
  default:
    std::cout << "WARNING: Unknown delta_pipi mode: " << delta_PiPi_mode << std::endl;
    std::cout << "Defaulting to delta_PiPi_DEFAULT." << std::endl;
    return delta_PiPi_DEFAULT(thermo, params);
    break;
  }
}

KOKKOS_INLINE_FUNCTION
double lambda_Pipi(const double* thermo, parameters params)
{
  int lambda_Pipi_mode = params.lambda_Pipi_mode;
  switch (lambda_Pipi_mode)
  {
  case LAMBDA_B_PIPI_DEFAULT:
    return lambda_Pipi_DEFAULT(thermo, params);
    break;
  case LAMBDA_B_PIPI_DISABLED:
    return 0.0; // No lambda PiPi contribution
    break;
  default:
    std::cout << "WARNING: Unknown lambda_Pipi mode: " << lambda_Pipi_mode << std::endl;
    std::cout << "Defaulting to lambda_Pipi_DEFAULT." << std::endl;
    return lambda_Pipi_DEFAULT(thermo, params);
    break;
  }
}

KOKKOS_INLINE_FUNCTION
double phi1(const double* thermo, parameters params)
{
  int phi1_mode = params.phi1_mode;
  switch (phi1_mode)
  {
  case PHI1_DEFAULT:
    return phi1_DEFAULT(thermo, params);
    break;
  case PHI1_DISABLED:
    return 0.0; // No phi1 contribution
    break;
  default:
    std::cout << "WARNING: Unknown phi1 mode: " << phi1_mode << std::endl;
    std::cout << "Defaulting to phi1_DEFAULT." << std::endl;
    return phi1_DEFAULT(thermo, params);
    break;
  }
}

double phi3(const double* thermo, parameters params)
{
  int phi3_mode = params.phi3_mode;
  switch (phi3_mode)
  {
  case PHI3_DEFAULT:
    return phi3_DEFAULT(thermo, params);
    break;
  case PHI3_DISABLED:
    return 0.0; // No phi3 contribution
    break;
  default:
    std::cout << "WARNING: Unknown phi3 mode: " << phi3_mode << std::endl;
    std::cout << "Defaulting to phi3_DEFAULT." << std::endl;
    return phi3_DEFAULT(thermo, params);
    break;
  }
}

double delta_pipi(const double* thermo, parameters params)
{
  int delta_pipi_mode = params.delta_pipi_mode;
  switch (delta_pipi_mode)
  {
  case DELTA_S_PIPI_DEFAULT:
    return delta_pipi_DEFAULT(thermo, params);
    break;
  case DELTA_S_PIPI_IS:
    return delta_pipi_IS(thermo, params);
    break;
  case DELTA_S_PIPI_DISABLED:
    return 0.0; // No delta PiPi contribution
    break;
  default:
    std::cout << "WARNING: Unknown delta_pipi mode: " << delta_pipi_mode << std::endl;
    std::cout << "Defaulting to delta_pipi_DEFAULT." << std::endl;
    return delta_pipi_DEFAULT(thermo, params);
    break;
  }
}

KOKKOS_INLINE_FUNCTION
double tau_pipi(const double* thermo, parameters params)
{
  int tau_pipi_mode = params.tau_pipi_mode;
  switch (tau_pipi_mode)
  {
  case TAU_PIPI_DEFAULT:
    return tau_pipi_DEFAULT(thermo, params);
    break;
  case TAU_PIPI_DISABLED:
    return 0.0; // No tau PiPi contribution
    break;
  default:
    std::cout << "WARNING: Unknown tau_pipi mode: " << tau_pipi_mode << std::endl;
    std::cout << "Defaulting to tau_pipi_DEFAULT." << std::endl;
    return tau_pipi_DEFAULT(thermo, params);
    break;
  }
}

KOKKOS_INLINE_FUNCTION
double lambda_piPi(const double* thermo, parameters params)
{
  int lambda_piPi_mode = params.lambda_piPi_mode;
  switch (lambda_piPi_mode)
  {
  case LAMBDA_S_PIPI_DEFAULT:
    return lambda_piPi_DEFAULT(thermo, params);
    break;
  case LAMBDA_S_PIPI_DISABLED:
    return 0.0; // No lambda PiPi contribution
    break;
  default:
    std::cout << "WARNING: Unknown lambda_piPi mode: " << lambda_piPi_mode << std::endl;
    std::cout << "Defaulting to lambda_piPi_DEFAULT." << std::endl;
    return lambda_piPi_DEFAULT(thermo, params);
    break;
  }
}

KOKKOS_INLINE_FUNCTION
double phi6(const double* thermo, parameters params)
{
  int phi6_mode = params.phi6_mode;
  switch (phi6_mode)
  {
  case PHI6_DEFAULT:
    return phi6_DEFAULT(thermo, params);
    break;
  case PHI6_DISABLED:
    return 0.0; // No phi6 contribution
    break;
  default:
    std::cout << "WARNING: Unknown phi6 mode: " << phi6_mode << std::endl;
    std::cout << "Defaulting to phi6_DEFAULT." << std::endl;
    return phi6_DEFAULT(thermo, params);
    break;
  }
}

KOKKOS_INLINE_FUNCTION
double phi7(const double* thermo, parameters params)
{
  int phi7_mode = params.phi7_mode;
  switch (phi7_mode)
  {
  case PHI7_DEFAULT:
    return phi7_DEFAULT(thermo, params);
    break;
  case PHI7_DISABLED:
    return 0.0; // No phi7 contribution
    break;
  default:
    std::cout << "WARNING: Unknown phi7 mode: " << phi7_mode << std::endl;
    std::cout << "Defaulting to phi7_DEFAULT." << std::endl;
    return phi7_DEFAULT(thermo, params);
    break;
  }
}





KOKKOS_INLINE_FUNCTION
Matrix<double, 3, 3> kappa(const double* thermo, parameters params)
{
  int diffusion_mode = params.diffusion_mode;
  switch (diffusion_mode)
  {
  case KAPPA_DEFAULT:
    return default_kappa(thermo, params);
    break;
  case KAPPA_DNMR:
    return Matrix<double, 3, 3>{0.0};
    std::cout << "DNMR kappa not implemented yet" << std::endl;
    break;
  case KAPPA_B_COTH:
    return kappaB_coth(thermo, params);
    break;
  case KAPPA_CONST_NB_OVER_T:
    return const_nB_over_T_kappa(thermo, params);
    break;
  case KAPPA_B_HOLO_TABLE:
    return kappaB_holo_table(thermo, params);
    break;
  default:
    return default_kappa(thermo, params);
    break;
  }
}

KOKKOS_INLINE_FUNCTION
Matrix<double, 3, 3> tauq_base(const double* thermo, parameters params);

/// tau_q of the configured relaxation mode, times the critical-slowing-down factor
/// (xi/xi_0)^p_tauq in xi_power mode (exactly 1 otherwise, so every other configuration is unchanged).
KOKKOS_INLINE_FUNCTION
Matrix<double, 3, 3> tauq(const double* thermo, parameters params)
{
  Matrix<double, 3, 3> tauq_matrix = tauq_base(thermo, params);
  const double slowing = critical_factor_tauq(thermo, params);
  if ( slowing != 1.0 )
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) tauq_matrix(i, j) *= slowing;
  return tauq_matrix;
}

KOKKOS_INLINE_FUNCTION
Matrix<double, 3, 3> tauq_base(const double* thermo, parameters params)
{

  if ( params.diffusion_relaxation_mode == TAU_Q_CATTANEO )
    return tauq_cattaneo(thermo, params);
  if ( params.diffusion_relaxation_mode == TAU_Q_KAPPA_OVER_N )
    return tauq_kappa_over_n(thermo, params);

  int diffusion_mode = params.diffusion_mode;
  switch (diffusion_mode)
  {
  case KAPPA_DEFAULT:
    return default_tauq(thermo, params);
    break;
  case KAPPA_DNMR:
    return Matrix<double, 3, 3>{0.0};
    std::cout << "DNMR tauq not implemented yet" << std::endl;
    break;
  case KAPPA_CONST_NB_OVER_T:
    return default_tauq(thermo, params);  // tau_q = 0.2/T
    break;
  case KAPPA_B_HOLO_TABLE:
    if (params.diffusion_relaxation_mode == TAU_Q_CB_OVER_T)
      return tauqB_only(thermo, params);   // tau_q = C_B/T
    else
      return default_tauq(thermo, params); // tau_q = 0.2/T
    break;
  case KAPPA_B_COTH:
    // tau_q = C_B/T (cb_over_T) or 0.2/T (constant_over_T, independent of C_B)
    if (params.diffusion_relaxation_mode == TAU_Q_CB_OVER_T)
      return tauqB_only(thermo, params);   // tau_q = C_B/T
    else
      return default_tauq(thermo, params); // tau_q = 0.2/T
    break;
  default:
    return default_tauq(thermo, params);
    break;
  }
}


//===============================
KOKKOS_INLINE_FUNCTION
double default_eta(const double *therm)
{
  double eta_over_s = 0.20;
  double w = therm[thermo_info::w];
  double T = therm[thermo_info::T];
  return eta_over_s*(w/T);
}

//===============================
KOKKOS_INLINE_FUNCTION
double constEta(const double *therm, const parameters params)
{
  double w = therm[thermo_info::w];
  double T = therm[thermo_info::T];
  if (T <= 0.0) return 0.0;
  return params.constant_eta_over_s*(w/T);
}

//===============================
KOKKOS_INLINE_FUNCTION
double JakiParam(const double *therm)
{
  // parameters hardcoded for now (TODO: make them settings)

  double T     = therm[thermo_info::T];
  double s     = therm[thermo_info::s];
  double TC    = 155.0/hbarc_MeVfm; // 173.9/197.3
  double z     = pow(0.66*T/TC,2);
  double alpha = 33./(12.*M_PI)*(z-1)/(z*log(z));
  return s*(0.0416762/pow(alpha,1.6)+ 0.0388977/pow(T,5.1) );
}

//===============================
KOKKOS_INLINE_FUNCTION
double LinearMusParam(const double *therm, const parameters params)
{
  // parameters hardcoded for now.. just to see how it works
  ///TODO: Create interface for changing parameters. For now, use constant eta/s as base viscosity
  double etaBase = params.constant_eta_over_s;
  double muSlope = 0.0033;
  double muB = therm[thermo_info::muB];
  double muS = therm[thermo_info::muS];
  double muQ = therm[thermo_info::muQ];
  double w = therm[thermo_info::w];
  double T = therm[thermo_info::T];
  if (T <= 0.0) return 0.0;
  return (etaBase + muSlope*(muB + muS + muQ))*(w/T);
}

//===============================
/// Piecewise-linear eta/s(muB) (BES Bayesian, PRC 110, 054905): nodes shear_eta_0/2/4
/// at muB = 0, 0.2, 0.4 GeV, constant beyond 0.4 GeV.
KOKKOS_INLINE_FUNCTION
double PiecewiseMuBParam(const double *therm, const parameters params)
{
  double w = therm[thermo_info::w];
  double T = therm[thermo_info::T];
  if (T <= 0.0) return 0.0;
  double muB_GeV = Kokkos::fabs(therm[thermo_info::muB])*hbarc_GeVfm;
  double eta_over_s;
  if (muB_GeV <= 0.2)
    eta_over_s = params.shear_eta_0
                 + (params.shear_eta_2 - params.shear_eta_0)*muB_GeV/0.2;
  else if (muB_GeV < 0.4)
    eta_over_s = params.shear_eta_2
                 + (params.shear_eta_4 - params.shear_eta_2)*(muB_GeV - 0.2)/0.2;
  else
    eta_over_s = params.shear_eta_4;
  return eta_over_s*(w/T);
}


//===============================
///TODO: add this in later..
KOKKOS_INLINE_FUNCTION
double InterpolantWrapper(const double *therm) { return 0.0; }


//==============================================================================
//==============================================================================
// Possible function choices for shear relaxation

//===============================
KOKKOS_INLINE_FUNCTION
double default_tau_pi(const double *therm, const parameters params) {
   return Kokkos::max( 5.0*eta(therm,params)/therm[thermo_info::w], 0.005 ); }

//===============================
KOKKOS_INLINE_FUNCTION
double tau_piGubser(const double *therm, const parameters params) {
   return (5.0*eta(therm,params))/therm[thermo_info::w];
   }

//===============================
KOKKOS_INLINE_FUNCTION
double tau_piMinval(const double *therm, const parameters params) {
   return Kokkos::max( (5.0*eta(therm,params))/therm[thermo_info::w], 0.001 ); }



//==============================================================================
KOKKOS_INLINE_FUNCTION
double delta_pipi_DEFAULT(const double *therm, const parameters params)
{
  // Default value for delta_PiPi

  return 4.*tau_pi(therm, params)/3.0;
}

KOKKOS_INLINE_FUNCTION
double delta_pipi_IS(const double *therm, const parameters params)
{
  return tau_pi(therm, params);
}


KOKKOS_INLINE_FUNCTION
double tau_pipi_DEFAULT(const double *therm, const parameters params)
{
  // Default value for tau_pipi
  return 10.*tau_pi(therm, params)/7.0;
}

KOKKOS_INLINE_FUNCTION
double lambda_piPi_DEFAULT(const double *therm, const parameters params)
{
  // Default value for lambda_piPi
  return 6./5.*tau_pi(therm, params);
}

double phi6_DEFAULT(const double *therm, const parameters params)
{
  // Default value for phi6
  return 0.;
}
KOKKOS_INLINE_FUNCTION
double phi7_DEFAULT(const double *therm, const parameters params)
{
  double pressure = therm[thermo_info::p];
  // Default value for phi7
  return 9./(70.*pressure);
}


//==============================================================================
//==============================================================================
// Possible function choices for zeta


//===============================
KOKKOS_INLINE_FUNCTION
double default_zeta(const double *therm)
{
  double zeta_over_s = 0.005;
  return zeta_over_s*(therm[thermo_info::w]/therm[thermo_info::T]);
}
//===============================
KOKKOS_INLINE_FUNCTION
double constZeta(const double *therm, const parameters params)
{
  double zeta_over_s = params.constant_zeta_over_s;
  double w = therm[thermo_info::w];
  double T = therm[thermo_info::T];
  return zeta_over_s*(w/T);
}

//===============================
KOKKOS_INLINE_FUNCTION
double zeta_DNMR_LeadingMass(const double *therm)
{
    //add this in later.. for now no bulk
    return 0.0;
}

//===============================
KOKKOS_INLINE_FUNCTION
double cs2_dependent_zeta(const double *therm, const parameters params)
{
  const double A = params.cs2_dependent_zeta_A;
  const double p = params.cs2_dependent_zeta_p;

  // optional tanh factor (with a power law) that suppresses the bulk at low T
  double factor = 1.0, th_x = 0.0, x_p = 0.0;
  if ( params.modulate_zeta_with_tanh )
  {
    const double T_transition = 150.0/constants::hbarc_MeVfm,
                 T_scale      = 10.0/constants::hbarc_MeVfm;
    th_x = tanh( ( therm[thermo_info::T] - T_transition ) / T_scale );
    x_p = pow(therm[thermo_info::T]/T_transition, p);
    factor = 0.5*(1.0 + x_p) + 0.5*(1.0 - x_p)*th_x;
  }

  const double zeta_over_s_local
                = A * factor * pow((1.0/3.0) - Kokkos::min(therm[thermo_info::cs2], 1.0), p);
  #ifdef DEBUG
  #ifndef __CUDACC__
   if ( therm[thermo_info::cs2] < 0.0 || therm[thermo_info::cs2] > 1.0 )
  {
    cout << "ERROR: " << zeta_over_s_local << "   "
        << x_p << "   "
        << th_x << "   "
        << factor << ";   "
        << therm[thermo_info::T] << "   "
        << therm[thermo_info::muB] << "   "
        << therm[thermo_info::muS] << "   "
        << therm[thermo_info::muQ] << "   "
        << therm[thermo_info::p] << "   "
        << therm[thermo_info::s] << "   "
        << therm[thermo_info::rhoB] << "   "
        << therm[thermo_info::rhoS] << "   "
        << therm[thermo_info::rhoQ] << "   "
        << therm[thermo_info::e] << "   "
        << therm[thermo_info::cs2] << "   "
        //<< therm.eos_name << "   "
        << A << "   "
        << p << "   "
        << pow((1.0/3.0)-therm[thermo_info::cs2], p) << endl;

    // do not tolerate this error
    if ( therm[thermo_info::cs2] < 0.0 )
    {
      cout << "cs2 went negative!" << endl;
      abort();
    }
  }
  if ( zeta_over_s_local > 0.5 )
  {
    cout << "LARGE ZETA/S: " << zeta_over_s_local << "   "
        << x_p << "   "
        << th_x << "   "
        << factor << ";   "
        << therm[thermo_info::T] << "   "
        << therm[thermo_info::muB] << "   "
        << therm[thermo_info::muS] << "   "
        << therm[thermo_info::muQ] << "   "
        << therm[thermo_info::p] << "   "
        << therm[thermo_info::s] << "   "
        << therm[thermo_info::rhoB] << "   "
        << therm[thermo_info::rhoS] << "   "
        << therm[thermo_info::rhoQ] << "   "
        << therm[thermo_info::e] << "   "
        << therm[thermo_info::cs2] << "   "
        //<< therm.eos_name << "   "
        << A << "   "
        << p << "   "
        << pow((1.0/3.0)-therm[thermo_info::cs2], p) << endl;
  }
  #endif
  #endif
  double w = therm[thermo_info::w];
  double T = therm[thermo_info::T];
  return zeta_over_s_local * w/T;
}

//==============================================================================
//==============================================================================
// Possible function choices for bulk relaxation

//===============================
// The tanh factor of cs2_dependent_zeta alone: all of zeta that survives in default_tau_Pi.
KOKKOS_INLINE_FUNCTION
double cs2_dependent_zeta_tanh_factor(const double *therm, const parameters params)
{
  if ( !params.modulate_zeta_with_tanh ) return 1.0;
  const double p = params.cs2_dependent_zeta_p;
  const double T_transition = 150.0/constants::hbarc_MeVfm,
               T_scale      = 10.0/constants::hbarc_MeVfm;
  const double th_x = tanh( ( therm[thermo_info::T] - T_transition ) / T_scale );
  const double x_p  = pow(therm[thermo_info::T]/T_transition, p);
  return 0.5*(1.0 + x_p) + 0.5*(1.0 - x_p)*th_x;
}

//===============================
KOKKOS_INLINE_FUNCTION
double default_tau_Pi(const double *therm, const parameters params)
{
  double p = params.cs2_dependent_zeta_p;
  const double causal_cs2 = therm[thermo_info::cs2] > 1 ? 1 : therm[thermo_info::cs2];

  if ( params.zeta_mode == ZETA_CS2_DEPENDENT )
  {
    // zeta = A f (1/3-cs2)^p w/T cancels: tau_Pi = coeff zeta/((1/3-cs2)^p w) = coeff A f/T
    double tau = params.bulk_tau_pi_coeff * params.cs2_dependent_zeta_A
                 * cs2_dependent_zeta_tanh_factor(therm, params)
                 / therm[thermo_info::T];
    if ( params.critical_scaling_bulk )
      tau *= critical_scaling_factor_bulk(therm, params);
    return Kokkos::max( tau, 0.1 );
  }

  // Every other zeta mode (constant / DNMR / default / interpolate) has no
  // (1/3-cs2)^p in zeta, so here the division -- and the guard -- are genuine.
  if ( pow((1.0/3.0-causal_cs2), p) < 1e-10 )
    return 1e10;
  else
    // tau_Pi = bulk_tau_pi_coeff zeta / ((1/3-cs2)^p w)
    return Kokkos::max( params.bulk_tau_pi_coeff*zeta(therm,params)/( pow((1.0/3.0-causal_cs2), p)*therm[thermo_info::w] ), 0.1 );
}

//===============================
// DNMR (14-moment) bulk relaxation time: tau_Pi = zeta / [15 (1/3 - cs2)^2 (e + P)]
KOKKOS_INLINE_FUNCTION
double tau_Pi_DNMR15(const double *therm, const parameters params)
{
  const double causal_cs2 = therm[thermo_info::cs2] > 1 ? 1 : therm[thermo_info::cs2];
  const double gap = Kokkos::max( 1.0/3.0 - causal_cs2, 1e-3 );
  const double w   = therm[thermo_info::w];
  if ( !(w > 0.0) ) return 0.1;
  const double tau = zeta(therm, params) / ( 15.0 * gap * gap * w );
  return Kokkos::max( tau, 0.1 );
}

KOKKOS_INLINE_FUNCTION
double delta_PiPi_DEFAULT(const double *therm, const parameters params)
{
  // Default value for delta_PiPi
  return 2.*tau_Pi(therm, params)/3.0;
}

KOKKOS_INLINE_FUNCTION
double delta_PiPi_IS(const double *therm, const parameters params)
{
  return tau_Pi(therm, params);
}

KOKKOS_INLINE_FUNCTION
double lambda_Pipi_DEFAULT(const double *therm, const parameters params)
{

  const double causal_cs2 = therm[thermo_info::cs2] > 1 ? 1 : therm[thermo_info::cs2];
  return 8.*(1./3. - causal_cs2)*tau_Pi(therm, params)/5.;
  
}

KOKKOS_INLINE_FUNCTION
double phi1_DEFAULT(const double *therm, const parameters params)
{
  // Default value for phi1
  return 0.0;
}
KOKKOS_INLINE_FUNCTION
double phi3_DEFAULT(const double *therm, const parameters params)
{
  // Default value for phi3
  return 0.0;
}


//===============================
// Diffusion parameters
KOKKOS_INLINE_FUNCTION
Matrix<double, 3, 3> default_kappa(const double *therm, const parameters params)
{
  Matrix<double, 3, 3> kappa_matrix;
  double T = therm[thermo_info::T];
  //construct sqrt(ni nj)/sqrt(mu_i mu_j) matrix (checking if the potential is not zero)
  double rhoQ = therm[thermo_info::rhoQ];
  double rhoS = therm[thermo_info::rhoS];
  double rhoB = therm[thermo_info::rhoB];
  Vector<double, 3> rho= {rhoB, rhoS, rhoQ};
  double muQ = therm[thermo_info::muQ];
  double muS = therm[thermo_info::muS];
  double muB = therm[thermo_info::muB];
  Vector<double, 3> mu= {muB, muS, muQ};

  Matrix<double, 3, 3> ni_nj;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      if (sqrt(abs(mu(i) * mu(j))) < TINY)
        ni_nj(i, j) = 0.0;
      else
        ni_nj(i, j) = sqrt(abs(rho(i) * rho(j))) / sqrt(abs(mu(i) * mu(j)));
  

  const double suppression = critical_factor_kappa(therm, params);

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      //kappa_matrix(i, j) = params.constant_kappa_over_T2[i][j]*ni_nj(i, j);
      kappa_matrix(i, j) = suppression * params.constant_kappa_over_T2[i][j]*(T*T);

  return kappa_matrix;
}


//===============================
// kappa_ij = const_ij sqrt(|rho_i rho_j|)/T (constant_kappa_over_T2 matrix): vanishes where
// the density does and has no coth(muB/T) edge amplification. tau_q = 0.2/T.
KOKKOS_INLINE_FUNCTION
Matrix<double, 3, 3> const_nB_over_T_kappa(const double *therm, const parameters params)
{
  Matrix<double, 3, 3> kappa_matrix{0.0};
  const double T = therm[thermo_info::T];
  if ( T <= TINY ) return kappa_matrix;
  Vector<double, 3> rho = { therm[thermo_info::rhoB],
                            therm[thermo_info::rhoS],
                            therm[thermo_info::rhoQ] };
  const double suppression = critical_factor_kappa(therm, params);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      kappa_matrix(i, j) = suppression * params.constant_kappa_over_T2[i][j]
                           * Kokkos::sqrt(Kokkos::abs(rho(i) * rho(j))) / T;
  return kappa_matrix;
}

//===============================
/// @brief Raises a prescribed tau_q (0.2/T or C_B/T) to the causality bound
/// tau_q >= D_B = kappa_BB/(T chi_BB) (Kapusta & Young, PRC 90, 044902, below Eq. 27).
KOKKOS_INLINE_FUNCTION
double tau_q_with_causality_floor(double tau_q, const double *therm,
                                  const parameters params)
{
  if ( !params.tau_q_causality_floor ) return tau_q;
  const double T      = therm[thermo_info::T];
  const double chi_BB = therm[thermo_info::chiBB];
  if ( !(T > TINY) || !(chi_BB > TINY) ) return tau_q;
  const double diffusion_constant = kappa(therm, params)(0, 0)/(T*chi_BB);
  if ( !Kokkos::isfinite(diffusion_constant) ) return tau_q;
  return Kokkos::fmax( tau_q, diffusion_constant );
};

//from PhysRevD.101.076007
KOKKOS_INLINE_FUNCTION
Matrix<double, 3, 3> default_tauq(const double *therm, const parameters params)
{
  Matrix<double, 3, 3> tauq_matrix;
  double T = therm[thermo_info::T];
  double ntot = therm[thermo_info::rhoB] + therm[thermo_info::rhoS] + therm[thermo_info::rhoQ];
  for (int i = 0; i < 3; ++i){
    for (int j = 0; j < 3; ++j) tauq_matrix(i, j) = 0.0;
    
  }

  //add diagonal terms
  const double tau_q = tau_q_with_causality_floor(0.2/T, therm, params);
  for (int i = 0; i < 3; ++i){
    tauq_matrix(i, i) += tau_q;
  }

  return tauq_matrix;
};

//===============================
// Baryon-only kappa_B (Eq. 16): kappa_B = (C_B/T) n_B [ (1/3) coth(mu_B/T) - n_B T/(e+P) ]
// (only the (B,B) entry; needs a B-only charge configuration).
KOKKOS_INLINE_FUNCTION
Matrix<double, 3, 3> kappaB_coth(const double *therm, const parameters params)
{
  Matrix<double, 3, 3> kappa_matrix{0.0};
  const double T   = therm[thermo_info::T];
  const double muB = therm[thermo_info::muB];
  const double nB  = therm[thermo_info::rhoB];
  const double w   = therm[thermo_info::w];   // enthalpy density = e + P

  // coth(muB/T) diverges as muB/T -> 0 (0*inf in near-empty cells): keep the argument off 0
  constexpr double muB_over_T_floor = 1.0e-2;
  const double x = ( T > TINY ) ? muB / T : 0.0;

  if ( T > TINY && Kokkos::abs(w) > TINY && Kokkos::abs(x) > muB_over_T_floor )
  {
    // |muB| and |nB|, so the coth term cannot flip sign at the table edge
    const double abs_x  = Kokkos::abs(x);
    const double abs_nB = Kokkos::abs(nB);
    // bracket > 0 at low density, < 0 in the dense core (diffusion suppressed there)
    const double bracket = (1.0/3.0)/Kokkos::tanh(abs_x) - abs_nB * T / w;
    // kappa_B = 0 once the density term wins (the core suppression)
    const double bracket_plus = Kokkos::fmax(bracket, 0.0);
    double kB = (params.C_B / T) * abs_nB * bracket_plus;
    // critical scaling (inert unless critical_scaling_diffusion is on)
    kB *= critical_factor_kappa(therm, params);
    // bracket_plus >= 0 already guarantees kB >= 0, so no hard clamp is needed -
    // just reject non-finite values.
    if ( Kokkos::isfinite(kB) )
      kappa_matrix(0, 0) = kB;
  }

  return kappa_matrix;
}

//===============================
// kappa_B = T sigma_B(T, mu_B) from the tabulated holographic EMD baryon conductivity
// (Rougemont, Noronha & Noronha-Hostler, arXiv:1507.06972), bilinear lookup. Outside the
// table (trustworthy for 130 <= T <= 300 MeV, mu_B <= 400 MeV) the edge value is used.
// Only the (B,B) entry; needs a B-only charge configuration.
KOKKOS_INLINE_FUNCTION
Matrix<double, 3, 3> kappaB_holo_table(const double *therm, const parameters params)
{
  Matrix<double, 3, 3> kappa_matrix{0.0};

  const int num_temperatures = params.holo_sigmaB_num_temperatures;
  const int num_potentials   = params.holo_sigmaB_num_potentials;
  if ( num_temperatures < 2 || num_potentials < 2 ) return kappa_matrix;

  const double T   = therm[thermo_info::T];
  const double muB = therm[thermo_info::muB];
  if ( !(T > TINY) ) return kappa_matrix;

  // fractional grid coordinates; the table is symmetric under muB -> -muB, so the
  // magnitude is what is looked up (matching kappaB_coth's symmetrized treatment).
  double temperature_coordinate = ( T - params.holo_sigmaB_temperature_min )
                                / params.holo_sigmaB_temperature_step;
  double potential_coordinate   = ( Kokkos::abs(muB) - params.holo_sigmaB_potential_min )
                                / params.holo_sigmaB_potential_step;

  temperature_coordinate = Kokkos::fmin( Kokkos::fmax(temperature_coordinate, 0.0),
                                         (double)(num_temperatures - 1) );
  potential_coordinate   = Kokkos::fmin( Kokkos::fmax(potential_coordinate, 0.0),
                                         (double)(num_potentials - 1) );

  const int iT  = Kokkos::min( (int)temperature_coordinate, num_temperatures - 2 );
  const int imu = Kokkos::min( (int)potential_coordinate,   num_potentials   - 2 );
  const double wT  = temperature_coordinate - (double)iT;
  const double wmu = potential_coordinate   - (double)imu;

  const double sigma_B =
      (1.0 - wT)*(1.0 - wmu)*params.holo_sigmaB(iT,     imu    )
    + (1.0 - wT)*       wmu *params.holo_sigmaB(iT,     imu + 1)
    +        wT *(1.0 - wmu)*params.holo_sigmaB(iT + 1, imu    )
    +        wT *       wmu *params.holo_sigmaB(iT + 1, imu + 1);

  // sigma_B is stored in natural units (fm^-1), so kappa_B comes out in fm^-2,
  // the same units as every other kappa mode.
  double kB = params.holo_diffusion_scale * T * sigma_B;

  // critical scaling: an extra factor on top of the EMD result (its CP lies outside the table)
  kB *= critical_factor_kappa(therm, params);

  if ( Kokkos::isfinite(kB) && kB > 0.0 )
    kappa_matrix(0, 0) = kB;

  return kappa_matrix;
}

//===============================
// Cattaneo baryon relaxation time (Kapusta & Young, PRC 90, 044902): tau_q = D_B / v_0^2 with
// D_B = kappa_BB/(T chi_BB) (= sigma_B/chi_BB for the holographic kappa_B) and v_0^2 the local
// cs2, or diffusion_wave_speed_sq if > 0. Clamped into [tau_q_min_over_T, tau_q_max_over_T]/T
// (cs2 -> 0 at a critical point, chi_BB -> 0 in the corona), then raised to the causality
// bound tau_q >= D_B. All three diagonals are filled so tau_q stays invertible.
KOKKOS_INLINE_FUNCTION
Matrix<double, 3, 3> tauq_cattaneo(const double *therm, const parameters params)
{
  Matrix<double, 3, 3> tauq_matrix{0.0};
  const double T = therm[thermo_info::T];
  if ( !(T > TINY) ) return tauq_matrix;

  const double tau_q_min = params.tau_q_min_over_T/T;
  const double tau_q_max = params.tau_q_max_over_T/T;

  // v_0^2: the configured value, or the local speed of sound squared. Anything
  // outside (0,1] would be acausal or nonsensical.
  double wave_speed_sq = params.diffusion_wave_speed_sq;
  if ( wave_speed_sq <= 0.0 ) wave_speed_sq = therm[thermo_info::cs2];

  const double chi_BB = therm[thermo_info::chiBB];
  double diffusion_constant = 0.0;
  double tau_q = tau_q_max;
  if ( wave_speed_sq > 0.0 && wave_speed_sq <= 1.0 && chi_BB > TINY )
  {
    // D_B = kappa_BB/(T chi_BB)  [Kapusta & Young Eq. (4): sigma = D dn/dmu]
    diffusion_constant = kappa(therm, params)(0, 0)/(T*chi_BB);
    tau_q = diffusion_constant/wave_speed_sq;   // tau = D/v_0^2
  }
  if ( !Kokkos::isfinite(tau_q) ) tau_q = tau_q_max;
  tau_q = Kokkos::fmin( Kokkos::fmax(tau_q, tau_q_min), tau_q_max );
  // Causality (Kapusta & Young, below their Eq. (27)): tau_q must exceed D_B, or
  // the large-k group velocity outruns light. This wins over the clamps above.
  if ( params.tau_q_causality_floor && Kokkos::isfinite(diffusion_constant) )
    tau_q = Kokkos::fmax( tau_q, diffusion_constant );

  for (int i = 0; i < 3; ++i) tauq_matrix(i, i) = tau_q;
  return tauq_matrix;
}

// tau_q = b_q kappa_BB / |n_B| (a time in fm)
KOKKOS_INLINE_FUNCTION
Matrix<double, 3, 3> tauq_kappa_over_n(const double *therm, const parameters params)
{
  Matrix<double, 3, 3> tauq_matrix{0.0};
  const double T = therm[thermo_info::T];
  if ( !(T > TINY) ) return tauq_matrix;

  const double tau_q_min = params.tau_q_min_over_T/T;
  const double tau_q_max = params.tau_q_max_over_T/T;

  const double kappa_BB = kappa(therm, params)(0, 0);
  const double n_B      = Kokkos::abs( therm[thermo_info::rhoB] );

  double tau_q = tau_q_max;
  if ( n_B > TINY && kappa_BB > 0.0 )
    tau_q = params.b_q*kappa_BB/n_B;
  if ( !Kokkos::isfinite(tau_q) ) tau_q = tau_q_max;
  tau_q = Kokkos::fmin( Kokkos::fmax(tau_q, tau_q_min), tau_q_max );

  // causality tau_q >= D_B, after the clamps (off with causality_floor: false)
  const double chi_BB = therm[thermo_info::chiBB];
  if ( params.tau_q_causality_floor && chi_BB > TINY )
  {
    const double diffusion_constant = kappa_BB/(T*chi_BB);
    if ( Kokkos::isfinite(diffusion_constant) )
      tau_q = Kokkos::fmax( tau_q, diffusion_constant );
  }

  for (int i = 0; i < 3; ++i) tauq_matrix(i, i) = tau_q;
  return tauq_matrix;
}

// tau_q = C_B / T on all three diagonals (only (B,B) is used; the rest keeps tau_q invertible).
KOKKOS_INLINE_FUNCTION
Matrix<double, 3, 3> tauqB_only(const double *therm, const parameters params)
{
  Matrix<double, 3, 3> tauq_matrix{0.0};
  const double T = therm[thermo_info::T];
  if ( T > TINY )
  {
    // C_B/T is a prescribed time, unconnected to kappa_B, so it can sit below the
    // Kapusta-Young bound - and does, with the holographic kappa_B. Impose it here.
    const double tau_q = tau_q_with_causality_floor(params.C_B/T, therm, params);
    for (int i = 0; i < 3; ++i)
      tauq_matrix(i, i) = tau_q;
  }
  return tauq_matrix;
};

//Critical scaling functions
KOKKOS_INLINE_FUNCTION
double critical_scaling_factor_bulk(const double *therm,const parameters params)
{
  double chiBB = therm[thermo_info::chiBB];
  double chiBB0 = therm[thermo_info::chiBB0];
  // chiBB < 0 in the unstable / gap region of a CP EoS (sqrt would be NaN): no enhancement
  if ( !(chiBB > 0.0) || !(chiBB0 > 0.0) ) return 1.0;
  double correlation_length      = Kokkos::sqrt(chiBB);
  double correlation_length_zero = Kokkos::sqrt(chiBB0);

  // xi_0: the table's local non-critical baseline chiBB_nc when present (ratio 1 off the
  // ridge), else the muB = 0 susceptibility with the global constant below.
  if ( !params.critical_scaling_local_xi0 )
    correlation_length_zero = 1.0176198446545879;

  if ( Kokkos::abs(correlation_length_zero - correlation_length) < TINY )
    return 1.0;

  double raw = 1.0 + pow((correlation_length/correlation_length_zero), 3);

  // local baseline: raw -> 2 off the ridge, so halve it to leave the bulk unchanged there
  if ( params.critical_scaling_local_xi0 ) raw *= 0.5;

  // The ridge is razor-thin: the raw factor reaches ~233 at the critical point and ~1259 at the
  // plane maximum, and it multiplies both zeta and tau_Pi.  Uncapped that is not survivable.
  if ( params.critical_scaling_max > 0.0 && raw > params.critical_scaling_max )
    raw = params.critical_scaling_max;

  return raw;
}

/// @brief Factor in [0, 1] suppressing diffusion near the critical point:
/// f = 1 - exp(-[(T-Tc)^2/(2 sigmaT^2) + (muB-muBc)^2/(2 sigmaB^2)]).
KOKKOS_INLINE_FUNCTION
double gaussian_suppression_diffusion(const double *therm, const parameters params)
{
  const double dT   = therm[thermo_info::T]   - params.critical_point_T;
  const double dmuB = therm[thermo_info::muB] - params.critical_point_muB;
  const double sigmaT   = params.critical_gaussian_width_T;
  const double sigmamuB = params.critical_gaussian_width_muB;
  const double exponent = 0.5 * ( dT*dT/(sigmaT*sigmaT) + dmuB*dmuB/(sigmamuB*sigmamuB) );
  return 1.0 - Kokkos::exp(-exponent);
}

/// @brief g(p) = min{[1 + (xi/xi_0)^p]/2, critical_scaling_max}, xi/xi_0 = sqrt(chi_BB/chi_BB,nc)
KOKKOS_INLINE_FUNCTION
double critical_xi_power_factor(const double *therm, const parameters params, const double power)
{
  const double chiBB  = therm[thermo_info::chiBB];
  const double chiBB0 = therm[thermo_info::chiBB0];
  if ( !(chiBB > 0.0) || !(chiBB0 > 0.0) ) return 1.0;
  double factor = 0.5*( 1.0 + pow( Kokkos::sqrt(chiBB/chiBB0), power ) );
  if ( !Kokkos::isfinite(factor) ) return 1.0;
  if ( params.critical_scaling_max > 0.0 && factor > params.critical_scaling_max )
    factor = params.critical_scaling_max;
  return factor;
}

/// @brief Critical factor on kappa_B: legacy Gaussian suppression, or the model-H
///        conductivity enhancement g(p_kappa) (xi_power, p = 1: kappa_n ~ xi). 1 when off.
KOKKOS_INLINE_FUNCTION
double critical_factor_kappa(const double *therm, const parameters params)
{
  if ( !params.critical_scaling_diffusion || !(therm[thermo_info::eos_type_is_table] > 0.5) )
    return 1.0;
  if ( params.critical_scaling_diffusion_mode == 1 )
    return critical_xi_power_factor(therm, params, params.critical_scaling_kappa_power);
  return gaussian_suppression_diffusion(therm, params);
}

/// @brief Critical slowing down tau_q = tau_q,nc g(p_tauq) (xi_power mode only; else 1).
KOKKOS_INLINE_FUNCTION
double critical_factor_tauq(const double *therm, const parameters params)
{
  if ( !params.critical_scaling_diffusion || params.critical_scaling_diffusion_mode != 1
    || !(therm[thermo_info::eos_type_is_table] > 0.5) )
    return 1.0;
  // the kappa-derived forms already carry g(p_kappa) (cattaneo also the (xi/xi_0)^2 of chi_BB):
  // take them back to tau_q,nc before applying g(p_tauq)
  double factor = critical_xi_power_factor(therm, params, params.critical_scaling_tauq_power);
  if ( params.diffusion_relaxation_mode == TAU_Q_CATTANEO
    || params.diffusion_relaxation_mode == TAU_Q_KAPPA_OVER_N )
  {
    const double g_kappa = critical_xi_power_factor(therm, params, params.critical_scaling_kappa_power);
    if ( g_kappa > 0.0 ) factor /= g_kappa;
    if ( params.diffusion_relaxation_mode == TAU_Q_CATTANEO )
    {
      const double chiBB  = therm[thermo_info::chiBB];
      const double chiBB0 = therm[thermo_info::chiBB0];
      if ( chiBB > 0.0 && chiBB0 > 0.0 && Kokkos::isfinite(chiBB/chiBB0) ) factor *= chiBB/chiBB0;
    }
  }
  return factor;
}

//==============================================================================
// Host-only loader for the holographic sigma_B table (arXiv:1507.06972 
// input/tables/holographic_sigmaB_EMD.dat): rows of T, mu_B, sigma_B [MeV] on a uniform
// grid; other lines are skipped. Converted to fm^-1 here.
inline void load_holographic_sigmaB_table(parameters & params,
                                          std::shared_ptr<Settings> settingsPtr)
{
  if ( params.diffusion_mode != KAPPA_B_HOLO_TABLE ) return;

  const std::string path = settingsPtr->holographic_sigmaB_table_path;
  std::ifstream table_file(path);
  if ( !table_file.is_open() )
  {
    std::cerr << "ERROR: could not open holographic sigma_B table '" << path
              << "'." << std::endl;
    Kokkos::finalize();
    exit(EXIT_FAILURE);
  }

  std::vector<double> temperatures, potentials, conductivities;
  std::string line;
  while ( std::getline(table_file, line) )
  {
    std::istringstream row(line);
    double T_MeV, muB_MeV, sigmaB_MeV;
    std::string leftover;
    if ( !(row >> T_MeV >> muB_MeV >> sigmaB_MeV) ) continue;  // header / blank
    if ( row >> leftover ) continue;                           // not a 3-column row
    temperatures.push_back(T_MeV);
    potentials.push_back(muB_MeV);
    conductivities.push_back(sigmaB_MeV);
  }
  table_file.close();

  if ( temperatures.empty() )
  {
    std::cerr << "ERROR: holographic sigma_B table '" << path
              << "' contains no 'T muB sigma_B' rows." << std::endl;
    Kokkos::finalize();
    exit(EXIT_FAILURE);
  }

  // Recover the two axes as the sorted unique values of each column.
  auto unique_sorted = [](std::vector<double> values)
  {
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    return values;
  };
  const std::vector<double> temperature_axis = unique_sorted(temperatures);
  const std::vector<double> potential_axis   = unique_sorted(potentials);

  const int num_temperatures = (int)temperature_axis.size();
  const int num_potentials   = (int)potential_axis.size();
  if ( num_temperatures < 2 || num_potentials < 2
       || (size_t)num_temperatures*(size_t)num_potentials != temperatures.size() )
  {
    std::cerr << "ERROR: holographic sigma_B table '" << path << "' is not a complete "
              << "rectangular grid (" << num_temperatures << " x " << num_potentials
              << " axes vs " << temperatures.size() << " rows)." << std::endl;
    Kokkos::finalize();
    exit(EXIT_FAILURE);
  }

  // The device lookup assumes uniform spacing on both axes; verify it here rather
  // than silently interpolating a table it cannot address.
  auto uniform_step = [&](const std::vector<double> & axis, const char * name)
  {
    const double step = (axis.back() - axis.front())/(double)(axis.size() - 1);
    for ( size_t i = 1; i < axis.size(); ++i )
      if ( std::abs( (axis[i] - axis[i-1]) - step ) > 1.0e-6*std::abs(step) )
      {
        std::cerr << "ERROR: holographic sigma_B table '" << path << "' has a "
                  << "non-uniform " << name << " axis at index " << i << "."
                  << std::endl;
        Kokkos::finalize();
        exit(EXIT_FAILURE);
      }
    return step;
  };
  const double temperature_step_MeV = uniform_step(temperature_axis, "temperature");
  const double potential_step_MeV   = uniform_step(potential_axis,   "chemical potential");

  // Fill the device view, converting MeV -> fm^-1 on every axis and on sigma_B.
  Kokkos::View<double**, MemorySpace> sigmaB_view("holographic_sigmaB",
                                                  num_temperatures, num_potentials);
  auto host_sigmaB = Kokkos::create_mirror_view(sigmaB_view);
  for ( size_t row = 0; row < temperatures.size(); ++row )
  {
    const int iT  = (int)std::lround( (temperatures[row] - temperature_axis.front())
                                      / temperature_step_MeV );
    const int imu = (int)std::lround( (potentials[row] - potential_axis.front())
                                      / potential_step_MeV );
    host_sigmaB(iT, imu) = conductivities[row]/hbarc_MeVfm;
  }
  Kokkos::deep_copy(sigmaB_view, host_sigmaB);

  params.holo_sigmaB                   = sigmaB_view;
  params.holo_sigmaB_num_temperatures  = num_temperatures;
  params.holo_sigmaB_num_potentials    = num_potentials;
  params.holo_sigmaB_temperature_min   = temperature_axis.front()/hbarc_MeVfm;
  params.holo_sigmaB_temperature_step  = temperature_step_MeV/hbarc_MeVfm;
  params.holo_sigmaB_potential_min     = potential_axis.front()/hbarc_MeVfm;
  params.holo_sigmaB_potential_step    = potential_step_MeV/hbarc_MeVfm;

  formatted_output::report("Loaded holographic sigma_B table '" + path + "': "
    + std::to_string(num_temperatures) + " x " + std::to_string(num_potentials)
    + " grid, T = " + std::to_string(temperature_axis.front()) + " - "
    + std::to_string(temperature_axis.back()) + " MeV, muB = "
    + std::to_string(potential_axis.front()) + " - "
    + std::to_string(potential_axis.back()) + " MeV; kappa_B = "
    + std::to_string(params.holo_diffusion_scale) + " * T * sigma_B");
}

}}

#endif