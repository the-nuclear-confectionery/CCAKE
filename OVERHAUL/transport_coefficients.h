#ifndef TRANSPORT_COEFFICIENTS_H
#define TRANSPORT_COEFFICIENTS_H
#include <string>
#include <vector>
#include <functional>
#include <math.h>
#include "settings.h"
#include "kernel.h"
#include "particle.h"
#include "new_matrix.h"
#include "eos.h"

using std::string;
using std::vector;

class TransportCoeficients
{

  public:
    TransportCoeficients();
    ~TransportCoeficients();

    double getEta();
    double getZeta();
    // Matrix getKappa();
    double getTauShear();
    double getTauBulk();
    // Matrix getTauDiffusive();
    void set_EquationOfStatePtr( EquationOfState * eosPtr_in );

    void initialize();


  private:
    string etaType, zetaType;
    string etaOption, zetaOption; // for the case of interp should we
    // read in path and load here or should I/O
    //load in directly?? Similar qeustion for EOS..

    double constEta();
    double eta_T_OV_w_IN;
    double JakiParam();
    double LinearMusParam();
    double InterpolantWrapper();
    function<double()> eta;
    
    double tauShearGubser();
    double tauShearMinval();
    function<double()> tauShear;

    double zeta_DNMR_LeadingMass();
    function<double()> zeta;

    double tauBulk_DNMR_LeadingMass();
    function<double()> tauBulk;


    EquationOfState * eosPtr = nullptr;

    









}



#endif
