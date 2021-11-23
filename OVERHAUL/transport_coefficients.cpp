#include <string>
#include <vector>
#include <functional>
#include <math.h>
#include "settings.h"
#include "kernel.h"
#include "particle.h"
#include "new_matrix.h"
#include "eos.h"
#include "transport_coefficients.h"

using std::string;
using std::vector;


//Constructors/destructor, may need to update later
TransportCoeficients::TransportCoeficients()
{

}
TransportCoeficients::~TransportCoeficients()
{

}

void TransportCoeficients::set_EquationOfStatePtr( EquationOfState * eosPtr_in )
{
  eosPtr = eosPtr_in;
}
///////////////////////////////////////////////////////////////
//Getter functions for eta, zeta, and their relaxation times//
/////////////////////////////////////////////////////////////
double TransportCoeficients::getEta()
{
    return eta();
}
double TransportCoeficients::getZeta()
{
    return zeta();
}
double TransportCoeficients::getTauShear()
{
    return tauShear();
}
double TransportCoeficients::getTauBulk()
{
    return tauBulk();
}

//////////////////////////////////////////////////////////////
//////////////possible function choices for eta//////////////
////////////////////////////////////////////////////////////
double TransportCoeficients::constEta()
{
    return eta_T_OV_w_IN*(eosPtr->w()/eosPtr->T())
}
double TransportCoeficients::JakiParam()
{
    // parameters hardcoded for now.. just to see how it works
    double TC=155; // 173.9/197.3
    double temp=eosPtr->T()*197.3/TC;
    double z=pow(0.66*temp,2);
    double alpha=33./(12.*PI)*(z-1)/(z*log(z));
    return eosPtr->s()*(0.0416762/pow(alpha,1.6)+ 0.0388977/pow(temp,5.1) );
}
double TransportCoeficients::LinearMusParam()
{
    // parameters hardcoded for now.. just to see how it works
    double etaBase = 0.08;
    double muSlope = 0.0033;
    return (etaBase + muSlope*(eosPtr->muB + 
    eosPtr->muS + eosPtr->muQ))*(eosPtr->w()/eosPtr->T())
}
double TransportCoeficients::InterpolantWrapper()
{
    //need interpolator for this, first
    return 0.0;
}

//////////////////////////////////////////////////////////////
////////possible function choices for shear relaxation///////
////////////////////////////////////////////////////////////
double TransportCoeficients::tauShearGubser()
{
    return (5*eta())/eosPtr->w();
}
double TransportCoeficients::tauShearMinval()
{
    double tau = (5*eta())/eosPtr->w();
    if (tau >= .001)
    {
        return tau
    }
    else
    {
        return 0.001
    }
}

//////////////////////////////////////////////////////////////
//////////////possible function choices for eta//////////////
////////////////////////////////////////////////////////////
double TransportCoeficients::zeta_DNMR_LeadingMass()
{
    //add this in later.. for now no bulk
    return 0.0;
}

//////////////////////////////////////////////////////////////
////////possible function choices for bulk relaxation///////
////////////////////////////////////////////////////////////
double TransportCoeficients::tauBulk_DNMR_LeadingMass()
{
    //add this in later.. for now no bulk
    return 0.0;
}
