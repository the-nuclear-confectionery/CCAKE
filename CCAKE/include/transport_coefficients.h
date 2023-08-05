#ifndef TRANSPORT_COEFFICIENTS_H
#define TRANSPORT_COEFFICIENTS_H

#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "eos.h"
#include "kernel.h"
#include "matrix.h"
#include "particle.h"
#include "settings.h"
#include "thermodynamic_info.h"

#include "Kokkos_Macros.hpp"

using std::string;
using std::vector;

namespace ccake{
class TransportCoefficients
{

  private:
    string etaMode, zetaMode;
    string shearRelaxMode, bulkRelaxMode;

    void initialize_eta(const string & etaMode_in);
    void initialize_tau_pi(const string & shearRelaxMode_in);
    void initialize_zeta(const string & zetaMode_in);
    void initialize_tau_Pi(const string & bulkRelaxMode_in);

    // Chris' defaults
    KOKKOS_INLINE_FUNCTION
    double default_eta(double *therm);
    KOKKOS_INLINE_FUNCTION
    double default_tau_pi(double *therm);
    KOKKOS_INLINE_FUNCTION
    double default_zeta(double *therm);
    KOKKOS_INLINE_FUNCTION
    double default_tau_Pi(double *therm);

    // other parmaterizations (not organized yet)
    KOKKOS_INLINE_FUNCTION
    double constEta(double *therm);
    double eta_T_OV_w_IN;
    KOKKOS_INLINE_FUNCTION
    double JakiParam(double *therm);
    KOKKOS_INLINE_FUNCTION
    double LinearMusParam(double *therm);
    KOKKOS_INLINE_FUNCTION
    double InterpolantWrapper(double *therm);
    
    KOKKOS_INLINE_FUNCTION
    double tau_piGubser(double *therm);
    KOKKOS_INLINE_FUNCTION
    double tau_piMinval(double *therm);

    KOKKOS_INLINE_FUNCTION
    double constZeta(double *therm);
    KOKKOS_INLINE_FUNCTION
    double zeta_DNMR_LeadingMass(double *therm);
    KOKKOS_INLINE_FUNCTION
    double cs2_dependent_zeta(double *therm);

    KOKKOS_INLINE_FUNCTION
    double tau_Pi_DNMR_LeadingMass(double *therm);

    std::shared_ptr<Settings> settingsPtr;


  public:
    TransportCoefficients() = delete; ///< Default constructor is deleted. Settings must be passed in.
    TransportCoefficients(std::shared_ptr<Settings> settingsPtr_in);
    ~TransportCoefficients(){}

    void set_SettingsPtr( std::shared_ptr<Settings> settingsPtr_in ) { settingsPtr = settingsPtr_in; }

    // these functions actually return the needed transport coefficients
    // It always receives an array with a set of particle properties and 
    // an integer with the number of entries in the array
    std::function<double(double*)> eta;
    std::function<double(double*)> tau_pi;
    std::function<double(double*)> zeta;
    std::function<double(double*)> tau_Pi;

};
}


#endif
