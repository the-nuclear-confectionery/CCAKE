#ifndef TRANSPORT_COEFFICIENTS_H
#define TRANSPORT_COEFFICIENTS_H

#include <cmath>
#include <string>

#include <Cabana_Core.hpp>

#include "eos.h"
#include "settings.h"
#include "thermodynamic_info.h"
#include "constants.h"

namespace ccake{
namespace transport_coefficients{
  enum {
    SHEAR_DEFAULT,
    SHEAR_CONST,
    SHEAR_JAKI,
    SHEAR_LINEAR_MUS,
    SHEAR_INTERPOLATOR
  };
  enum {
    TAU_PI_SHEAR_DEFAULT,
    TAU_PI_SHEAR_MINVAL,
    TAU_PI_SHEAR_GUBSER,
    TAU_PI_SHEAR_CONSTANT
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

  struct parameters
  {
    int shear_mode, zeta_mode, shear_relaxation_mode, bulk_relaxation_mode;
    double constant_eta_over_s, constant_zeta_over_s;
    double constant_tau_pi;
    double cs2_dependent_zeta_A, cs2_dependent_zeta_p;
    bool modulate_zeta_with_tanh;
  };

  KOKKOS_INLINE_FUNCTION
  parameters setup_parameters(std::shared_ptr<Settings> settingsPtr);
  KOKKOS_INLINE_FUNCTION
  double eta(const double* thermo, parameters params);
  KOKKOS_INLINE_FUNCTION
  double tau_pi(const double* thermo, parameters params);
  KOKKOS_INLINE_FUNCTION
  double zeta(const double* thermo, parameters params);
  KOKKOS_INLINE_FUNCTION
  double tau_Pi(const double* thermo, parameters params);

  // Chris' defaults
  KOKKOS_INLINE_FUNCTION
  double default_eta(const double *therm);
  KOKKOS_INLINE_FUNCTION
  double default_tau_pi(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  double default_zeta(const double *therm);
  KOKKOS_INLINE_FUNCTION
  double default_tau_Pi(const double *therm, const parameters params);

  // other parmaterizations (not organized yet)
  KOKKOS_INLINE_FUNCTION
  double constEta(const double *therm, const parameters params);
  KOKKOS_INLINE_FUNCTION
  double JakiParam(const double *therm);
  KOKKOS_INLINE_FUNCTION
  double LinearMusParam(const double *therm, const parameters params);
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
  double tau_Pi_DNMR_LeadingMass(const double *therm);

}}


#endif
