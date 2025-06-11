#ifndef TRANSPORT_COEFFICIENTS_CPP
#define TRANSPORT_COEFFICIENTS_CPP
#include "transport_coefficients.h"

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
  else
  {
    std::cout << "WARNING: Unknown diffusion mode: " << diffusionMode << std::endl;
    std::cout << "Defaulting to constant_over_T2." << std::endl;
    params.diffusion_mode = KAPPA_DEFAULT; // default to KAPPA_DEFAULT if not specified
  }

  std::string diffusionRelaxMode = settingsPtr->diffusionRelaxMode;
  if (diffusionRelaxMode == "constant_over_T")
  {
    params.diffusion_relaxation_mode = TAU_Q_DEFAULT;
  }
  else{ 
    std::cout << "WARNING: Unknown diffusion relaxation mode: " << diffusionRelaxMode << std::endl
              << "Defaulting to constant_over_T." << std::endl;
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
  params.constant_zeta_over_s = settingsPtr->constant_zeta_over_s;
  params.cs2_dependent_zeta_A = settingsPtr->cs2_dependent_zeta_A;
  params.cs2_dependent_zeta_p = settingsPtr->cs2_dependent_zeta_p;
  params.modulate_zeta_with_tanh = settingsPtr->modulate_zeta_with_tanh;
  params.constant_kappa_over_T2 = settingsPtr->constant_kappa_over_T2;

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
  switch (zeta_mode)
  {
  case ZETA_DEFAULT:
    return default_zeta(thermo);
    break;
  case ZETA_CONSTANT:
    return constZeta(thermo, params);
    break;
  case ZETA_DNMR:
    return zeta_DNMR_LeadingMass(thermo);
    break;
  case ZETA_CS2_DEPENDENT:
    return cs2_dependent_zeta(thermo, params);
    break;
  default:
    return default_zeta(thermo);
    break;
  }
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
    return tau_Pi_DNMR_LeadingMass(thermo);
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
  default:
    return default_kappa(thermo, params);
    break;
  }
}

KOKKOS_INLINE_FUNCTION
Matrix<double, 3, 3> tauq(const double* thermo, parameters params)
{
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
  return therm[thermo_info::s] * eta_over_s;
}

//===============================
KOKKOS_INLINE_FUNCTION
double constEta(const double *therm, const parameters params)
{
  double w = therm[thermo_info::w];
  double T = therm[thermo_info::T];
  return params.constant_eta_over_s*(w/T);
}

//===============================
KOKKOS_INLINE_FUNCTION
double JakiParam(const double *therm)
{
  //picked the easiest one with functional dependence
  //parameters hardcoded for now. just to see how it works
  ///TODO: Create interface for changing parameters

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

  return (etaBase + muSlope*(muB + muS + muQ))*(w/T);
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
  return 6./5.;
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
  return therm[thermo_info::s] * zeta_over_s;
}

//===============================
KOKKOS_INLINE_FUNCTION
double constZeta(const double *therm, const parameters params)
{
  double zeta_over_s = params.constant_zeta_over_s;
  return therm[thermo_info::s] * zeta_over_s;
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

  //----------------------------------------------
  //!!!!!  Bulk is too large in low-T regime
  //!!!!!  ==>> add modulating tanh-factor with
  //!!!!!       power-law dependence to suppress
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
  if ( zeta_over_s_local > 0.1 )
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

  return zeta_over_s_local*therm[thermo_info::s];
}

//==============================================================================
//==============================================================================
// Possible function choices for bulk relaxation

//===============================
KOKKOS_INLINE_FUNCTION
double default_tau_Pi(const double *therm, const parameters params)
{
  double p = params.cs2_dependent_zeta_p;
  const double causal_cs2 = therm[thermo_info::cs2] > 1 ? 1 : therm[thermo_info::cs2];
  if ( pow((1.0/3.0-causal_cs2), p) < 1e-10 )
    return 1e10;
  else
    return Kokkos::max( 5.0*zeta(therm,params)/( pow((1.0/3.0-causal_cs2), p)*therm[thermo_info::w] ), 0.1 );
}

//===============================
KOKKOS_INLINE_FUNCTION
double tau_Pi_DNMR_LeadingMass(const double *therm) { return 0.0; }

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
  

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      //kappa_matrix(i, j) = params.constant_kappa_over_T2[i][j]*ni_nj(i, j);
      kappa_matrix(i, j) = params.constant_kappa_over_T2[i][j]*(T*T);
      
  return kappa_matrix;
}


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
  for (int i = 0; i < 3; ++i){
    tauq_matrix(i, i) += 0.2/T;
  }
      
  return tauq_matrix;
};

}}
#endif