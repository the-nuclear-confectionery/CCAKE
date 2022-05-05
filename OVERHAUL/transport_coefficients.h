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

  public:
    TransportCoefficients();
    ~TransportCoefficients();

    void set_SettingsPtr( Settings * settingsPtr_in );

    void initialize(  const string & etaType_in      = "default",
                      const string & tauShearType_in = "default",
                      const string & zetaType_in     = "default",
                      const string & tauBulkType_in  = "default" );

    void setTherm(thermodynamic_info & thermo_from_particle);

    double getEta();
    double getZeta();
    double getTauShear();
    double getTauBulk();


  private:
    thermodynamic_info therm;
    string etaType, zetaType;
    string etaOption, zetaOption; /*for the case of interp should we
    read in path and load here or should I/O
    load in directly?? Similar qeustion for EOS..*/
    string tauShearType, tauBulkType;

    // Chris' defaults
    double default_eta();
    double default_tauShear();
    double default_zeta();
    double default_tauBulk();

    // other parmaterizations (not organized yet)
    double constEta();
    double eta_T_OV_w_IN;
    double JakiParam();
    double LinearMusParam();
    double InterpolantWrapper();
    double NoShear();
    function<double()> eta;
    
    double tauShearGubser();
    double tauShearMinval();
    function<double()> tauShear;

    double zeta_DNMR_LeadingMass();
    double NoBulk();
    function<double()> zeta;

    double tauBulk_DNMR_LeadingMass();
    function<double()> tauBulk;

    Settings * settingsPtr   = nullptr;

};



#endif
