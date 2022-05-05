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

    void initialize();

    void setTherm(thermodynamic_info & thermo_from_particle); //maybe constructor??
    // change input to unordered map???? I think better that way so trcoef
    // doesnt need to be friends with particle

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

    static double constEta();
    static double eta_T_OV_w_IN;
    static double JakiParam();
    static double LinearMusParam();
    static double InterpolantWrapper();
    static double NoShear();
    function<double()> eta;
    
    static double tauShearGubser();
    static double tauShearMinval();
    function<double()> tauShear;

    static double zeta_DNMR_LeadingMass();
    static double NoBulk();
    function<double()> zeta;

    static double tauBulk_DNMR_LeadingMass();
    function<double()> tauBulk;

    Settings * settingsPtr   = nullptr;

};



#endif
