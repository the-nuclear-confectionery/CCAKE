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
  if (diffusionMode == "default")
  {
    params.diffusion_mode = KAPPA_DEFAULT;
  }
  else if (diffusionMode == "DNMR")
  {
    params.diffusion_mode = KAPPA_DNMR;
  }

  params.constant_eta_over_s = settingsPtr->constant_eta_over_s;
  params.constant_zeta_over_s = settingsPtr->constant_zeta_over_s;
  params.cs2_dependent_zeta_A = settingsPtr->cs2_dependent_zeta_A;
  params.cs2_dependent_zeta_p = settingsPtr->cs2_dependent_zeta_p;
  params.modulate_zeta_with_tanh = settingsPtr->modulate_zeta_with_tanh;
  params.kappa_matrix = settingsPtr->kappa_matrix;

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
      //kappa_matrix(i, j) = params.kappa_matrix[i][j]*ni_nj(i, j);
      kappa_matrix(i, j) = params.kappa_matrix[i][j]*(T*T)/0.716;
      
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