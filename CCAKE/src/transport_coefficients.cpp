#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "constants.h"
#include "eos.h"
#include "kernel.h"
#include "matrix.h"
#include "particle.h"
#include "settings.h"
#include "transport_coefficients.h"

using std::string;
using std::vector;
using namespace ccake;

/// @brief Constructor for TransportCoefficients class
/// @param etaType_in String with description of shear viscosity to be used
/// @param tau_piType_in String with description of shear relaxation time to be used
/// @param zetaType_in String with description of bulk viscosity to be used
/// @param tau_PiType_in String with description of bulk relaxation time to be used
TransportCoefficients::TransportCoefficients(std::shared_ptr<Settings> settingsPtr_in):
settingsPtr(settingsPtr_in)
{
  etaMode = settingsPtr->etaMode;
  shearRelaxMode = settingsPtr->shearRelaxMode;
  zetaMode = settingsPtr->zetaMode;
  bulkRelaxMode = settingsPtr->bulkRelaxMode;
  
  // Set shear viscosity
  initialize_eta( etaMode );

  // Set bulk viscosity
  initialize_zeta( zetaMode );

  // Set shear relaxation
  initialize_tau_pi( shearRelaxMode );

  // Set bulk relaxation
  initialize_tau_Pi( bulkRelaxMode );
}


//==============================================================================
void TransportCoefficients::initialize_eta(const string & etaType_in)
{
  // set chosen eta parameterization
  etaMode = etaType_in;

  // assign corresponding function
  if (etaMode == "default")
  {
    eta = [this](const double *therm){ return default_eta(therm); };
  }
  else if (etaMode == "constant")
  {
    eta_T_OV_w_IN = settingsPtr->constant_eta_over_s;
    eta = [this](const double *therm){ return constEta(therm); };
  }
  else if (etaMode == "JakiParam")
  {
    eta = [this](const double *therm){ return JakiParam(therm); };
  }
  else if (etaMode == "LinearMus")
  {
    eta = [this](const double *therm){ return LinearMusParam(therm); };
  }
  else if (etaMode == "interpolate")
  {
    //use parameter to find directory of table, then
    // execute interpolation
    eta = [this](const double *therm){ return InterpolantWrapper(therm); };
  }
  else
  {
    cout << "Shear viscosity specification " << etaMode << " not recognized. Now exiting.\n";
    exit(1);
  }
}


//==============================================================================
void TransportCoefficients::initialize_tau_pi(const string & tau_piType_in)
{
  shearRelaxMode = tau_piType_in;

  if (shearRelaxMode == "default")
  {
    tau_pi = [this](const double *therm){ return default_tau_pi(therm); };
  }
  else if (shearRelaxMode == "minVal")
  {
    tau_pi = [this](const double *therm){ return tau_piMinval(therm); };
  }
  else if (shearRelaxMode == "Gubser")
  {
    /* these consistency checks maybe should be done in settings.h? */
    if (etaMode != "constant")
    {
      std::cout << "Shear viscosity must be constant for Gubser. "
              "Check Input_Parameters.  Now exiting.\n";
      exit(1);
    }
    if (zetaMode != "constant" || abs(settingsPtr->constant_zeta_over_s) > 1e-6)
    {
      std::cout << "You have chosen zetaMode = " << zetaMode
                << " and zeta/s = " << abs(settingsPtr->constant_zeta_over_s) 
                << ".\n" << "Bulk viscosity must be zero"
                   " for Gubser. Check Input_Parameters.  "
                   " Now exiting.\n";
      exit(1);
    }
    tau_pi = [this](const double *therm){ return tau_piGubser(therm); };
  }
  else
  {
    cout << "Tau shear specification " << shearRelaxMode << " not recognized. Now exiting.\n";
    exit(1);
  }
}


//==============================================================================
void TransportCoefficients::initialize_zeta(const string & zetaType_in)
{
  zetaMode     = zetaType_in;

  if (zetaMode == "default")
  {
    zeta = [this](const double *therm){ return default_zeta(therm); };
  }
  else if (zetaMode == "constant")
  {
    zeta = [this](const double *therm){ return constZeta(therm); };
  }
  else if (zetaMode == "DNMR")
  {
    zeta = [this](const double *therm){ return zeta_DNMR_LeadingMass(therm); };
  }
  else if (zetaMode == "interpolate")
  {
    //use parameter to find directory of table, then
    // execute interpolation
    zeta = [this](const double *therm){ return InterpolantWrapper(therm); };
  }
  else if (zetaMode == "cs2_dependent")
  {
    zeta = [this](const double *therm){ return cs2_dependent_zeta(therm); };
  }
  else
  {
    cout << "Bulk viscosity specification " << zetaMode << " not recognized. Now exiting.\n";
    exit(1);
  }
}


//==============================================================================
void TransportCoefficients::initialize_tau_Pi(const string & tau_PiType_in)
{
  bulkRelaxMode  = tau_PiType_in;

  if (bulkRelaxMode == "default")
  {
    tau_Pi = [this](const double *therm){ return default_tau_Pi(therm); };
  }
  else if (bulkRelaxMode == "DNMR")
  {
    tau_Pi = [this](const double *therm){ return tau_Pi_DNMR_LeadingMass(therm); };
  }
  else 
  {
    cout << "Tau bulk specification " << bulkRelaxMode << " not recognized. Now exiting.\n";
    exit(1);   
  }
}

//===============================
double TransportCoefficients::default_eta(const double *therm)
{
  double eta_over_s = 0.20;
  return therm[thermo_info::s] * eta_over_s;
}

//===============================
double TransportCoefficients::constEta(const double *therm)
{
  double w = therm[thermo_info::w];
  double T = therm[thermo_info::T];
  return eta_T_OV_w_IN*(w/T);
}

//===============================
double TransportCoefficients::JakiParam(const double *therm)
{
  //picked the easiest one with functional dependence
  // parameters hardcoded for now.. just to see how it works

  double T     = therm[thermo_info::T];
  double s     = therm[thermo_info::s];
  double TC    = 155.0/hbarc_MeVfm; // 173.9/197.3
  double z     = pow(0.66*T/TC,2);
  double alpha = 33./(12.*PI)*(z-1)/(z*log(z));
  return s*(0.0416762/pow(alpha,1.6)+ 0.0388977/pow(T,5.1) );
}

//===============================
double TransportCoefficients::LinearMusParam(const double *therm)
{
  // parameters hardcoded for now.. just to see how it works
  double etaBase = 0.08;
  double muSlope = 0.0033;
  double muB = therm[thermo_info::muB];
  double muS = therm[thermo_info::muS];
  double muQ = therm[thermo_info::muQ];
  double w = therm[thermo_info::w];
  double T = therm[thermo_info::T];

  return (etaBase + muSlope*(muB + muS + muQ))*(w/T);
}

//===============================
///TODO: add this in later..
double TransportCoefficients::InterpolantWrapper(const double *therm) { return 0.0; }


//==============================================================================
//==============================================================================
// Possible function choices for shear relaxation

//===============================
double TransportCoefficients::default_tau_pi(const double *therm) {
   return std::max( 5.0*eta(therm)/therm[thermo_info::w], 0.005 ); }

//===============================
double TransportCoefficients::tau_piGubser(const double *therm) {
   return (5.0*eta(therm))/therm[thermo_info::w];
   }

//===============================
double TransportCoefficients::tau_piMinval(const double *therm) {
   return std::max( (5.0*eta(therm))/therm[thermo_info::w], 0.001 ); }


//==============================================================================
//==============================================================================
// Possible function choices for zeta


//===============================
double TransportCoefficients::default_zeta(const double *therm)
{
  double zeta_over_s = 0.005;
  return therm[thermo_info::s] * zeta_over_s;
}

//===============================
double TransportCoefficients::constZeta(const double *therm)
{
  double zeta_over_s = settingsPtr->constant_zeta_over_s;
  return therm[thermo_info::s] * zeta_over_s;
}

//===============================
double TransportCoefficients::zeta_DNMR_LeadingMass(const double *therm)
{
    //add this in later.. for now no bulk
    return 0.0;
}


//===============================
double TransportCoefficients::cs2_dependent_zeta(const double *therm)
{
  const double A = settingsPtr->cs2_dependent_zeta_A;
  const double p = settingsPtr->cs2_dependent_zeta_p;

  //----------------------------------------------
  //!!!!!  Bulk is too large in low-T regime
  //!!!!!  ==>> add modulating tanh-factor with
  //!!!!!       power-law dependence to suppress
  double factor = 1.0, th_x = 0.0, x_p = 0.0;
  if ( settingsPtr->modulate_zeta_with_tanh )
  {
    const double T_transition = 150.0/constants::hbarc_MeVfm,
                 T_scale      = 10.0/constants::hbarc_MeVfm;
    th_x = tanh( ( therm[thermo_info::T] - T_transition ) / T_scale );
    x_p = pow(therm[thermo_info::T]/T_transition, p);
    factor = 0.5*(1.0 + x_p) + 0.5*(1.0 - x_p)*th_x;
  }

  const double zeta_over_s_local
                = A * factor * pow((1.0/3.0) - std::min(therm[thermo_info::cs2], 1.0), p);
//  cout << "Check zeta/s: "
//        << therm.T*hbarc_MeVfm << "   "
//        << zeta_over_s_local << endl;

  if ( therm[thermo_info::cs2] < 0.0 || therm[thermo_info::cs2] > 1.0 )
  {
    cout << "ERROR: " << zeta_over_s_local << "   "
        << x_p << "   "
        << th_x << "   "
        << factor << ";   "
        << therm[thermo_info::T] << "   "
        << therm[thermo_info::muB] << "   "
        << therm[thermo_info::muS] << "   "
        << therm[thermo_info::muQ] << "   "
        << therm[thermo_info::p] << "   "
        << therm[thermo_info::s] << "   "
        << therm[thermo_info::rhoB] << "   "
        << therm[thermo_info::rhoS] << "   "
        << therm[thermo_info::rhoQ] << "   "
        << therm[thermo_info::e] << "   "
        << therm[thermo_info::cs2] << "   "
        //<< therm.eos_name << "   "
        << A << "   "
        << p << "   "
        << pow((1.0/3.0)-therm[thermo_info::cs2], p) << endl;

    // do not tolerate this error
    if ( therm[thermo_info::cs2] < 0.0 )
    {
      cout << "cs2 went negative!" << endl;
      abort();
    }
  }
  if ( zeta_over_s_local > 0.1 )
  {
    cout << "LARGE ZETA/S: " << zeta_over_s_local << "   "
        << x_p << "   "
        << th_x << "   "
        << factor << ";   "
        << therm[thermo_info::T] << "   "
        << therm[thermo_info::muB] << "   "
        << therm[thermo_info::muS] << "   "
        << therm[thermo_info::muQ] << "   "
        << therm[thermo_info::p] << "   "
        << therm[thermo_info::s] << "   "
        << therm[thermo_info::rhoB] << "   "
        << therm[thermo_info::rhoS] << "   "
        << therm[thermo_info::rhoQ] << "   "
        << therm[thermo_info::e] << "   "
        << therm[thermo_info::cs2] << "   "
        //<< therm.eos_name << "   "
        << A << "   "
        << p << "   "
        << pow((1.0/3.0)-therm[thermo_info::cs2], p) << endl;
  }

  return zeta_over_s_local*therm[thermo_info::s];
}

//==============================================================================
//==============================================================================
// Possible function choices for bulk relaxation

//===============================
double TransportCoefficients::default_tau_Pi(const double *therm)
{
//cout << "inside check: " << therm.cs2 << "   " << zeta() << "   " << therm.w << endl;
  if ( (1.0/3.0-therm[thermo_info::cs2])*(1.0/3.0-therm[thermo_info::cs2]) < 1e-10 )
    return 1e10;
  else
    return std::max( 5.0*zeta(therm)/(pow((1.0/3.0-therm[thermo_info::cs2]),2.0)*therm[thermo_info::w]), 0.1 );
}

//===============================
double TransportCoefficients::tau_Pi_DNMR_LeadingMass(const double *therm) { return 0.0; }
