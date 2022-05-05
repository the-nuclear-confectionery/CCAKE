#include <string>
#include <vector>
#include <functional>
#include <math.h>
#include "settings.h"
#include "kernel.h"
#include "particle.h"
#include "matrix.h"
#include "eos.h"
#include "transport_coefficients.h"

using std::string;
using std::vector;


//Constructors/destructor, may need to update later
TransportCoefficients::TransportCoefficients()
{
    etaType = settingsPtr->eta;
    etaOption = settingsPtr->etaOption;
    zetaType = settingsPtr->zeta;
    zetaOption = settingsPtr->zetaOption;
}
TransportCoefficients::~TransportCoefficients()
{

}

// void TransportCoefficients::set_EquationOfStatePtr( EquationOfState * eosPtr_in )
// {
//   eosPtr = eosPtr_in;
// }
void TransportCoefficients::set_SettingsPtr( Settings * settingsPtr_in )
{
  settingsPtr = settingsPtr_in;
}

//Setter for thermodynamic information
void setTherm(thermodynamic_info & thermo_from_particle)
{
    therm = thermo_from_particle;
}

///////////////////////////////////////////////////////////////
//Getter functions for eta, zeta, and their relaxation times//
/////////////////////////////////////////////////////////////
double TransportCoefficients::getEta()
{
    return eta();
}
double TransportCoefficients::getZeta()
{
    return zeta();
}
double TransportCoefficients::getTauShear()
{
    return tauShear();
}
double TransportCoefficients::getTauBulk()
{
    return tauBulk();
}

//////////////////////////////////////////////////////////////
//////////////possible function choices for eta//////////////
////////////////////////////////////////////////////////////
double TransportCoefficients::constEta()
{
    double w = therm.w;
    double T = therm.T;
    return eta_T_OV_w_IN*(w/T);
}
double TransportCoefficients::JakiParam()
{
    //picked the easiest one with functional dependence
    // parameters hardcoded for now.. just to see how it works

    double T = therm.T;
    double s = therm.s;
    double TC=155; // 173.9/197.3
    double TovTC = (T*hbarc_MeVfm)/TC;
    double z=pow(0.66*TovTC,2);
    double alpha=33./(12.*PI)*(z-1)/(z*log(z));
    return s*(0.0416762/pow(alpha,1.6)+ 0.0388977/pow(temp,5.1) );
}
double TransportCoefficients::LinearMusParam()
{
    // parameters hardcoded for now.. just to see how it works
    double etaBase = 0.08;
    double muSlope = 0.0033;
    double muB = therm.muB;
    double muS = therm.muS;
    double muQ = therm.muQ;
    double w = therm.w;
    double T = therm.T;

    return (etaBase + muSlope*(muB + muS + muQ))*(w/T);
}
double TransportCoefficients::InterpolantWrapper()
{
    //need interpolator for this, first
    return 0.0;
}
double TransportCoefficients::NoShear()
{
    return 0.0; // no shear means no shear
}

//////////////////////////////////////////////////////////////
////////possible function choices for shear relaxation///////
////////////////////////////////////////////////////////////
double TransportCoefficients::tauShearGubser()
{
    w = therm.w;
    return (5*eta(therm))/w;
}
double TransportCoefficients::tauShearMinval()
{
    w = therm.w;
    double tau = (5*eta(therm))/w;
    if (tau >= .001)
    {
        return tau;
    }
    else
    {
        return 0.001;
    }
}

//////////////////////////////////////////////////////////////
//////////////possible function choices for zeta//////////////
////////////////////////////////////////////////////////////
double TransportCoefficients::zeta_DNMR_LeadingMass()
{
    //add this in later.. for now no bulk
    return 0.0;
}
double TransportCoefficients::NoBulk()
{
    return 0.0; // of course, conformal zeta returns 0
}

//////////////////////////////////////////////////////////////
////////possible function choices for bulk relaxation///////
////////////////////////////////////////////////////////////
double TransportCoefficients::tauBulk_DNMR_LeadingMass()
{
    //add this in later.. for now no bulk
    return 0.0;
}

//////////////////////////////////////////////////
//////INITIALIZE THE TRANSPORT COEFFICIENTS//////
////////////////////////////////////////////////
void TransportCoefficients::initialize()
{
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
    else if (etaType == "LinearMus")
    {
        eta = LinearMusParam;
    }
    else if (etaType = "interpolate")
    {
        //use etaOption to find directory of table, then
        // execute interpolation
        eta = InterpolantWrapper;
    }
    else if (etaType = "NoShear")
    {
        eta = NoShear;
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
        /* these consistency checks maybe should be done in settings.h? */
        if (etaType != "constant")
        {
            cout << "Shear viscosity must be constant 
            for Gubser. Check Input_Parameters.  
            Now exiting." << endl;
            exit(1);
        }
        if (zetaType != "NoBulk")
        {
            cout << "Bulk viscosity must be conformal 
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
    else if (zetaType = "NoBulk")
    {
        zeta = NoBulk;
    }
    else if (zetaType = "interpolate")
    {
        //use zetaOption to find directory of table, then
        // execute interpolation
        zeta = InterpolantWrapper;
    }
    else
    {
        cout << "Bulk viscosity specification not recognized. Now exiting." << endl;
        exit(1);
    }

// SET BULK RELAXATION
    if (tauBulkType == "Default")
    {
        tauBulk = tauBulk_DNMR_LeadingMass;
    }
    else 
    {
        cout << "Tau bulk specification not 
        recognized. Now exiting." << endl;
        exit(1);   
    }


}