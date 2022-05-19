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

using std::string;
using std::vector;

class TransportCoefficients
{

  private:
    thermodynamic_info therm;
    string etaType, zetaType;
    string tau_piType, tau_PiType;

    void initialize_eta(const string & etaType_in);
    void initialize_tau_pi(const string & tau_piType_in);
    void initialize_zeta(const string & zetaType_in);
    void initialize_tau_Pi(const string & tau_PiType_in);

    // Chris' defaults
    double default_eta();
    double default_tau_pi();
    double default_zeta();
    double default_tau_Pi();

    // other parmaterizations (not organized yet)
    double constEta();
    double eta_T_OV_w_IN;
    double JakiParam();
    double LinearMusParam();
    double InterpolantWrapper();
    double NoShear();
    
    double tau_piGubser();
    double tau_piMinval();

    double constZeta();
    double zeta_DNMR_LeadingMass();
    double NoBulk();
    double cs2_dependent_zeta();

    double tau_Pi_DNMR_LeadingMass();

    Settings * settingsPtr   = nullptr;


  public:
    TransportCoefficients(){}
    ~TransportCoefficients(){}

    void set_SettingsPtr( Settings * settingsPtr_in ) { settingsPtr = settingsPtr_in; }

    // use this to set combinations of TC parameterizations automatically
    void initialize( const string & mode = "default" );
    void initialize( const string & etaType_in,  const string & tau_piType_in,
                     const string & zetaType_in, const string & tau_PiType_in );

    void setTherm( thermodynamic_info & thermo_from_particle ) { therm = thermo_from_particle; }

    // these functions actually return the needed transport coefficients
    std::function<double()> eta;
    std::function<double()> tau_pi;
    std::function<double()> zeta;
    std::function<double()> tau_Pi;

};



#endif
