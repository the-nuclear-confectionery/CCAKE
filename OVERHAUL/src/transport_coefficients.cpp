#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "../include/constants.h"
#include "../include/eos.h"
#include "../include/kernel.h"
#include "../include/matrix.h"
#include "../include/particle.h"
#include "../include/settings.h"
#include "../include/transport_coefficients.h"

using std::string;
using std::vector;

//==============================================================================
//==============================================================================
// INITIALIZE THE TRANSPORT COEFFICIENTS

//===============================
// batch initialization
void TransportCoefficients::initialize( const string & mode )
{
  if ( mode == "default" )
    initialize( "default", "default", "default", "default" );
  else if ( mode == "Gubser" )
    initialize( "constant", "Gubser", "NoBulk", "default" );
  else
  {
    std::cout << "TransportCoefficients::mode value = "
              << mode << " not supported!\n";
    exit(8);
  }
}


//===============================
// explicit (individual) initialization
void TransportCoefficients::initialize( const string & etaType_in,
                                        const string & tau_piType_in,
                                        const string & zetaType_in,
                                        const string & tau_PiType_in )
{
  cout << "Using etaType_in = " << etaType_in << endl;
  cout << "Using tau_piType_in = " << tau_piType_in << endl;
  cout << "Using zetaType_in = " << zetaType_in << endl;
  cout << "Using tau_PiType_in = " << tau_PiType_in << endl;

  // Set shear viscosity
  initialize_eta( etaType_in );

  // Set shear relaxation
  initialize_tau_pi( tau_piType_in );

  // Set bulk viscosity
  initialize_zeta( zetaType_in );

  // Set bulk relaxation
  initialize_tau_Pi( tau_PiType_in );
}


//==============================================================================
void TransportCoefficients::initialize_eta(const string & etaType_in)
{
  // set chosen eta parameterization
  etaType = etaType_in;

  // assign corresponding function
  if (etaType == "default")
  {
    eta = [this]{ return default_eta(); };
  }
  else if (etaType == "constant")
  {
    eta_T_OV_w_IN = settingsPtr->constant_eta_over_s;
    eta = [this]{ return constEta(); };
  }
  else if (etaType == "JakiParam")
  {
    eta = [this]{ return JakiParam(); };
  }
  else if (etaType == "LinearMus")
  {
    eta = [this]{ return LinearMusParam(); };
  }
  else if (etaType == "interpolate")
  {
    //use parameter to find directory of table, then
    // execute interpolation
    eta = [this]{ return InterpolantWrapper(); };
  }
  else if (etaType == "NoShear")
  {
    eta = [this]{ return NoShear(); };
  }
  else
  {
    cout << "Shear viscosity specification " << etaType << " not recognized. Now exiting.\n";
    exit(1);
  }
}


//==============================================================================
void TransportCoefficients::initialize_tau_pi(const string & tau_piType_in)
{
  tau_piType = tau_piType_in;

  if (tau_piType == "default")
  {
    tau_pi = [this]{ return default_tau_pi(); };
  }
  else if (tau_piType == "minVal")
  {
    tau_pi = [this]{ return tau_piMinval(); };
  }
  else if (tau_piType == "Gubser")
  {
    /* these consistency checks maybe should be done in settings.h? */
    if (etaType != "constant")
    {
      std::cout << "Shear viscosity must be constant for Gubser. "
              "Check Input_Parameters.  Now exiting.\n";
      exit(1);
    }
    if (zetaType != "NoBulk")
    {
      std::cout << "Bulk viscosity must be zero"
                   " for Gubser. Check Input_Parameters.  "
                   " Now exiting.\n";
      exit(1);
    }
    tau_pi = [this]{ return tau_piGubser(); };
  }
  else
  {
    cout << "Tau shear specification " << tau_piType << " not recognized. Now exiting.\n";
    exit(1);
  }
}


//==============================================================================
void TransportCoefficients::initialize_zeta(const string & zetaType_in)
{
  zetaType     = zetaType_in;

  if (zetaType == "default")
  {
    zeta = [this]{ return default_zeta(); };
  }
  else if (zetaType == "constant")
  {
    zeta = [this]{ return constZeta(); };
  }
  else if (zetaType == "DNMR")
  {
    zeta = [this]{ return zeta_DNMR_LeadingMass(); };
  }
  else if (zetaType == "NoBulk")
  {
    zeta = [this]{ return NoBulk(); };
  }
  else if (zetaType == "interpolate")
  {
    //use parameter to find directory of table, then
    // execute interpolation
    zeta = [this]{ return InterpolantWrapper(); };
  }
  else if (zetaType == "cs2_dependent")
  {
    zeta = [this]{ return cs2_dependent_zeta(); };
  }
  else
  {
    cout << "Bulk viscosity specification " << zetaType << " not recognized. Now exiting.\n";
    exit(1);
  }
}


//==============================================================================
void TransportCoefficients::initialize_tau_Pi(const string & tau_PiType_in)
{
  tau_PiType  = tau_PiType_in;

  if (tau_PiType == "default")
  {
    tau_Pi = [this]{ return default_tau_Pi(); };
  }
  else if (tau_PiType == "DNMR")
  {
    tau_Pi = [this]{ return tau_Pi_DNMR_LeadingMass(); };
  }
  else 
  {
    cout << "Tau bulk specification " << tau_PiType << " not recognized. Now exiting.\n";
    exit(1);   
  }
}




//==============================================================================
//==============================================================================
// Possible function choices for eta

//===============================
double TransportCoefficients::default_eta()
{
  double eta_over_s = 0.20;
  return therm.s * eta_over_s;
}

//===============================
double TransportCoefficients::constEta()
{
  double w = therm.w;
  double T = therm.T;
  return eta_T_OV_w_IN*(w/T);
}

//===============================
double TransportCoefficients::JakiParam()
{
  //picked the easiest one with functional dependence
  // parameters hardcoded for now.. just to see how it works

  double T     = therm.T;
  double s     = therm.s;
  double TC    = 155.0/hbarc_MeVfm; // 173.9/197.3
  double z     = pow(0.66*T/TC,2);
  double alpha = 33./(12.*PI)*(z-1)/(z*log(z));
  return s*(0.0416762/pow(alpha,1.6)+ 0.0388977/pow(T,5.1) );
}

//===============================
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

//===============================
double TransportCoefficients::InterpolantWrapper() { return 0.0; }

//===============================
double TransportCoefficients::NoShear() { return 0.0; }


//==============================================================================
//==============================================================================
// Possible function choices for shear relaxation

//===============================
double TransportCoefficients::default_tau_pi() { return std::max( 5.0*eta()/therm.w, 0.005 ); }

//===============================
double TransportCoefficients::tau_piGubser() { return (5.0*eta())/therm.w; }

//===============================
double TransportCoefficients::tau_piMinval() { return std::max( (5.0*eta())/therm.w, 0.001 ); }



//==============================================================================
//==============================================================================
// Possible function choices for zeta


//===============================
double TransportCoefficients::default_zeta()
{
  double zeta_over_s = 0.005;
  return therm.s * zeta_over_s;
}

//===============================
double TransportCoefficients::constZeta()
{
  double zeta_over_s = settingsPtr->constant_zeta_over_s;
  return therm.s * zeta_over_s;
}

//===============================
double TransportCoefficients::zeta_DNMR_LeadingMass()
{
    //add this in later.. for now no bulk
    return 0.0;
}

//===============================
double TransportCoefficients::NoBulk() { return 0.0; }

//===============================
double TransportCoefficients::cs2_dependent_zeta()
{
//  const double A = 8.0*constants::pi/15.0;
//  const double p = 2.0;
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
    th_x = tanh( ( therm.T - T_transition ) / T_scale );
    x_p = pow(therm.T/T_transition, p);
    factor = 0.5*(1.0 + x_p) + 0.5*(1.0 - x_p)*th_x;
  }

  cout << "Check zeta/s: " << therm.T*hbarc_MeVfm << "   "
        << A*factor*pow((1.0/3.0)-therm.cs2, p) << endl;
  const double zeta_over_s_local = A*factor*pow((1.0/3.0)-therm.cs2, p);

  if ( zeta_over_s_local > 100.0 )
  {
    cout << zeta_over_s_local << "   "
        << x_p << "   "
        << th_x << "   "
        << factor << ";   "
        << therm.T << "   "
        << therm.muB << "   "
        << therm.muS << "   "
        << therm.muQ << "   "
        << therm.p << "   "
        << therm.s << "   "
        << therm.rhoB << "   "
        << therm.rhoS << "   "
        << therm.rhoQ << "   "
        << therm.e << "   "
        << therm.cs2 << "   "
        << therm.eos_name << "   "
        << A << "   "
        << p << "   "
        << pow((1.0/3.0)-therm.cs2, p) << endl;
    abort();
}

  return A*factor*pow((1.0/3.0)-therm.cs2, p)*therm.s;
}

//==============================================================================
//==============================================================================
// Possible function choices for bulk relaxation

//===============================
double TransportCoefficients::default_tau_Pi()
{
//cout << "inside check: " << therm.cs2 << "   " << zeta() << "   " << therm.w << endl;
  if ( (1.0/3.0-therm.cs2)*(1.0/3.0-therm.cs2) < 1e-10 )
    return 1e10;
  else
    return std::max( 5.0*zeta()/(pow((1.0/3.0-therm.cs2),2.0)*therm.w), 0.1 );
}

//===============================
double TransportCoefficients::tau_Pi_DNMR_LeadingMass() { return 0.0; }
