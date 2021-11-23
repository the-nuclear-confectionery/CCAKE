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

//////////////////////////////////////////////////
//////INITIALIZE THE TRANSPORT COEFFICIENTS//////
////////////////////////////////////////////////
void TransportCoefficients::initialize()
{
    // assuming I/O has already set the necessary
    // parameters

    //SET SHEAR VISCOSITY
    if (etaType == "constant")
    {
        eta_T_OV_w_IN = stod(etaOption);
        eta = constEta;
    }
    else if (etaType == "JakiParam")
    {
        eta = JakiParam;
    }
    else if (etaType == "linearMus")
    {
        eta = LinearMusParam;
    }
    else if (etaType = "interpolate")
    {
        //use etaOption to find directory of table, then
        // execute interpolation
        eta = InterpolantWrapper;
    }
    else
    {
        cout << "Shear viscosity specification not 
        recognized. Now exiting." << endl;
        exit(1);
    }

// SET SHEAR RELAXATION
    if (tauShearType == "Default")
    {
        tauShear = tauShearMinval;
    }
    else if (tauShearType == "Gubser")
    {
        if (etaType != "constant")
        {
            cout << "Shear viscosity must be constant 
            for Gubser. Check Input_Parameters.  
            Now exiting." << endl;
            exit(1);
        }
        tauShear = tauShearGubser;
    }
    else
    {
        cout << "Tau shear specification not 
        recognized. Now exiting." << endl;
        exit(1);
    }

// SET BULK VISCOSITY
    if (zetaType = "Default")
    {
        zeta = zeta_DNMR_LeadingMass;
    }
    else if (zetaType = "interpolate")
    {
        //use zetaOption to find directory of table, then
        // execute interpolation
        zeta = InterpolantWrapper;
    }
    else
    {
        cout << "Bulk viscosity specification not 
        recognized. Now exiting." << endl;
        exit(1);
    }

// SET BULK RELAXATION
    if (tauBulkType == "Default")
    {
        zeta = tauBulk_DNMR_LeadingMass;
    }
    else 
    {
        cout << "Tau bulk specification not 
        recognized. Now exiting." << endl;
        exit(1);   
    }


}