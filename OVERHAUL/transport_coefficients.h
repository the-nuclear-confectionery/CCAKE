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

    void setTherm(thermodynamic_info & thermo_from_particle); //maybe constructor??
    // change input to unordered map???? I think better that way so trcoef
    // doesnt need to be friends with particle

    double getEta();// delete for others
    double getZeta();
    // Matrix getKappa();
    double getTauShear();
    double getTauBulk();
    // Matrix getTauDiffusive();
    // void set_EquationOfStatePtr( EquationOfState * eosPtr_in );
    void set_SettingsPtr( Settings * settingsPtr_in );

    void initialize();


  private:
    thermodynamic_info therm;
    string etaType, zetaType;
    string etaOption, zetaOption; /*for the case of interp should we
    read in path and load here or should I/O
    load in directly?? Similar qeustion for EOS..*/
    string tauShearType, tauBulkType;

    double constEta();
    double eta_T_OV_w_IN;
    double JakiParam();
    double LinearMusParam();
    double InterpolantWrapper();
    double NoShear();
    function<double(thermodynamic_info & therm)> eta;
    
    double tauShearGubser();
    double tauShearMinval();
    function<double()> tauShear;

    double zeta_DNMR_LeadingMass();
    double NoBulk();
    function<double()> zeta;

    double tauBulk_DNMR_LeadingMass();
    function<double()> tauBulk;


    // EquationOfState * eosPtr = nullptr;
    Settings * settingsPtr   = nullptr;


};



#endif
