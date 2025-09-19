#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "../include/constants.h"
#include "../include/eos.h"
#include "../include/eos_base.h"
#include "../include/eos_conformal.h"
#include "../include/eos_conformal_diagonal.h"
#include "../include/eos_header.h"
#include "../include/eos_nonconformal_extension.h"
#include "../include/eos_table.h"
#include "../include/eos_tanh_conformal.h"
#include "../include/formatted_output.h"
#include "../include/rootfinder.h"
#include "eos_conformal_cosh.h"

using namespace constants;
using namespace ccake;

using std::vector;
using std::string;
using std::to_string;

////////////////////////////////////////////////////////////////////////////////
void EquationOfState::set_SettingsPtr( std::shared_ptr<Settings> settingsPtr_in ) { settingsPtr = settingsPtr_in; }

////////////////////////////////////////////////////////////////////////////////
/// @brief Initializes the equation of state.
/// @details This function initializes the equation of state by setting up the
/// chosen EOSs and reading in the corresponding equation of state tables if
/// needed. It also includes an option to run a closure test on the equation
/// of state.
/// TODO: allow user to request closure test from input file
/// @return void
/// @see EquationOfState::set_up_chosen_EOSs()
/// @see EquationOfState::run_closure_test()
void EquationOfState::init()
{
  formatted_output::report("Initializing equation of state");

  formatted_output::update("Reading in equation of state tables");
  eos_path = settingsPtr->eos_path; ///\TODO: Is it really necessary?

  set_up_chosen_EOSs();

  #ifdef DEBUG
  bool do_eos_checks = settingsPtr->perform_eos_checks;
  if ( do_eos_checks )
    run_closure_test();
  #endif
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Initializes the conformal equation of state.
/// @details This function initializes the conformal equation of state by
/// setting up the corresponding EOS object and adding it to the chosen_EOSs
/// vector. The expression for the conformal equation of state is given by
/// \f{equation}{
///    p(T,\mu_B, \mu_Q, \mu_S)
///  = c T_0^4 \left[ \left( \frac{T}{T_0} \right)^2
///       +\sum_{X=B,S,Q} \left( \frac{\mu_X}{\mu_{X,0}} \right)^2 \right]^2
/// \f}
/// @param c The constant c in the conformal equation of state.
/// @param T0 The reference temperature \f$T_0\f$ in the conformal equation of state.
/// @param muB0 The reference baryon chemical potential \f$\mu_{B,0}\f$ in the conformal
/// equation of state.
/// @param muQ0 The reference electric charge chemical potential \f$\mu_{Q,0}\f$ in the
/// conformal equation of state.
/// @param muS0 The reference strangeness chemical potential \f$\mu_{S,0}\f$ in the
/// conformal equation of state.
/// @return void
/// @see EquationOfState::set_up_chosen_EOSs()
void EquationOfState::init_conformal(const double c, const double T0,
                                     const double muB0, const double muQ0, 
                                     const double muS0){
    formatted_output::update("Setting up conformal equation of state");
    //formatted_output::comment(
    //  "This conformal equation of state treats all T and mu axes equivalently "
    //  "and assumes an ideal gas of massless gluons, \"2.5\" massless quark "
    //  "flavors, and Nc = 3 colors.  Quadratic cross-terms are included.");

    // set minima and maxima (can be arbitrarily large)
    vector<double> tbqs_minima = { 0.0,          -TBQS_INFINITY, -TBQS_INFINITY, -TBQS_INFINITY };
    vector<double> tbqs_maxima = { TBQS_INFINITY, TBQS_INFINITY,  TBQS_INFINITY,  TBQS_INFINITY };

    // add EoS to vector
    chosen_EOSs.push_back( std::make_shared<EoS_conformal>(
                            c, T0, muB0, muS0, muQ0,
                            tbqs_minima, tbqs_maxima ) );

    formatted_output::detail("Set up conforma with following parameters:");
    formatted_output::detail( "c    = " + to_string(c) );
    formatted_output::detail( "T0   = " + to_string(T0) );
    formatted_output::detail( "muB0 = " + to_string(muB0) );
    formatted_output::detail( "muQ0 = " + to_string(muQ0) );
    formatted_output::detail( "muS0 = " + to_string(muS0) );

    default_eos_name = "conformal";
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Initializes the conformal diagonal equation of state.
/// @details This function initializes the conformal diagonal equation of state
/// by setting up the corresponding EOS object and adding it to the chosen_EOSs
/// vector. The expression for the conformal diagonal equation of state is given
/// by \f{equation}{
///    p(T, \mu_B, \mu_Q, \mu_S)
///     = c T_0^4 \left[ \left( \frac{T}{T_0} \right)^4
///          +\sum_{X=B,S,Q} \left( \frac{\mu_X}{\mu_{X,0}} \right)^4\right] ,
/// \f}
/// @param c The constant c in the conformal equation of state.
/// @param T0 The reference temperature \f$T_0\f$ in the conformal equation of state.
/// @param muB0 The reference baryon chemical potential \f$\mu_{B,0}\f$ in the conformal
/// equation of state.
/// @param muQ0 The reference electric charge chemical potential \f$\mu_{Q,0}\f$ in the
/// conformal equation of state.
/// @param muS0 The reference strangeness chemical potential \f$\mu_{S,0}\f$ in the
/// conformal equation of state.
/// @return void
/// @see EquationOfState::set_up_chosen_EOSs()
void EquationOfState::init_conformal_diagonal(const double c, const double T0,
                             const double muB0, const double muQ0, 
                             const double muS0){
    formatted_output::update("Setting diagonal conformal equation of state");
    //formatted_output::comment(
    //  "This conformal equation of state treats all T and mu axes equivalently "
    //  "and assumes an ideal gas of massless gluons, \"2.5\" massless quark "
    //  "flavors, and Nc = 3 colors.  Only quartic (diagonal) terms are "
    //  "included.");

    // set minima and maxima for rootfinder (can be arbitrarily large)
    vector<double> tbqs_minima = { 0.0,          -TBQS_INFINITY, -TBQS_INFINITY, -TBQS_INFINITY };
    vector<double> tbqs_maxima = { TBQS_INFINITY, TBQS_INFINITY,  TBQS_INFINITY,  TBQS_INFINITY };

    // add matched conformal EoS to vector of EoSs
    chosen_EOSs.push_back( std::make_shared<EoS_conformal_diagonal>(
                            c, T0, muB0, muS0, muQ0,
                            tbqs_minima, tbqs_maxima, "conformal_diagonal" ) );
    
    formatted_output::detail("Set up conformal diagonal with following parameters:");
    formatted_output::detail( "c    = " + to_string(c) );
    formatted_output::detail( "T0   = " + to_string(T0) );
    formatted_output::detail( "muB0 = " + to_string(muB0) );
    formatted_output::detail( "muQ0 = " + to_string(muQ0) );
    formatted_output::detail( "muS0 = " + to_string(muS0) );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Initializes the conformal diagonal equation of state.
/// @details This function initializes the tanh modulated "conformal" equation 
/// of state by setting up the corresponding EOS object and adding it to the 
/// chosen_EOSs vector. The expression for this equation of state is given
/// by \f{equation}{
/// p(T,\mu_B, \mu_Q, \mu_S)
///         = c T_0^4 \tanh\left(\frac{T-T_c}{T_s}\right) 
///         \times \left[ \left( \frac{T}{T_0} \right)^2
///          +\sum_{X=B,S,Q} \left( \frac{\mu_X}{\mu_{X,0}} \right)^2 \right]^2
/// \f}
/// @param c The constant c in the conformal equation of state.
/// @param T0 The reference temperature \f$T_0\f$ in the conformal equation of state.
/// @param muB0 The reference baryon chemical potential \f$\mu_{B,0}\f$ in the conformal
/// equation of state.
/// @param muQ0 The reference electric charge chemical potential \f$\mu_{Q,0}\f$ in the
/// conformal equation of state.
/// @param muS0 The reference strangeness chemical potential \f$\mu_{S,0}\f$ in the
/// conformal equation of state.
/// @param Tc Control the position of the modulating tanh function
/// @param Ts Control the strength of the modulating tanh function
/// @return void
/// @see EquationOfState::set_up_chosen_EOSs()
void EquationOfState::init_tanh_conformal(const double c, const double T0,
                         const double muB0, const double muQ0, 
                         const double muS0, const double Tc, const double Ts)
{
      formatted_output::update("setting tanh-modulated \"conformal\" equation "
                               "of state as fallback");
      formatted_output::detail("all coefficients matched to p/T^4 at grid limits");

      // set minima and maxima for rootfinder (can be arbitrarily large)
      vector<double> tbqs_minima = { 0.0,          -TBQS_INFINITY, -TBQS_INFINITY, -TBQS_INFINITY };
      vector<double> tbqs_maxima = { TBQS_INFINITY, TBQS_INFINITY,  TBQS_INFINITY,  TBQS_INFINITY };

      formatted_output::detail( "set up tanh-modulated \"conformal\" with following parameters:");
      formatted_output::detail( "c    = " + to_string(c) );
      formatted_output::detail( "T0   = " + to_string(T0) );
      formatted_output::detail( "muB0 = " + to_string(muB0) );
      formatted_output::detail( "muQ0 = " + to_string(muQ0) );
      formatted_output::detail( "muS0 = " + to_string(muS0) );
      formatted_output::detail( "Tc   = " + to_string(Tc) );
      formatted_output::detail( "Ts   = " + to_string(Ts) );

      // add matched conformal EoS to vector of EoSs
      chosen_EOSs.push_back( std::make_shared<EoS_tanh_conformal>(
                              c, T0, muB0, muS0, muQ0, Tc, Ts,
                              tbqs_minima, tbqs_maxima, "tanh_conformal" ) );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Sets up the equation of state based on the settings provided.
/// @details This function initializes the equation of state based on the
/// eos_type provided in the settings file. It supports three types of equation
/// of state: conformal, conformal diagonal, and table.
///
/// The chosen_EOSs vector is populated with the corresponding EOS object based
/// on the eos_type.
void EquationOfState::set_up_chosen_EOSs()
{
	tbqsPosition.resize(4);

  //============================================================================
  // THIS IF-CHAIN INITIALIZES THE DEFAULT EOS TO USE
  //============================================================================
  // SET UP CONFORMAL EOS
  if ( settingsPtr->eos_type == "conformal" )
  {
    const double Nc = 3.0, Nf = 2.5;  // u+d massless, s 'half massless'
    double c  = pi*pi*(2.0*(Nc*Nc-1.0)+(7.0/2.0)*Nc*Nf)/90.0;
    double T0 = 1.0, muB0 = 1.0, muQ0 = 1.0, muS0 = 1.0; // trivial scales
    init_conformal(c, T0, muB0, muQ0, muS0);
    init_conformal_diagonal(c, T0, muB0, muQ0, muS0);
    default_eos_name = "conformal";
  } else if ( settingsPtr->eos_type == "cosh_conformal" )
  {
    formatted_output::update("Setting up conformal cosh equation of state");
    //formatted_output::comment(
    //  "This conformal equation of state treats all T and mu axes equivalently "
    //  "and assumes an ideal gas of massless gluons, \"2.5\" massless quark "
    //  "flavors, and Nc = 3 colors.  Quadratic cross-terms are included.");

    double c  = 1/M_PI/M_PI;
    double T0 = 1.0, muB0 = 1.0, muQ0 = 1.0, muS0 = 1.0; // trivial scales

    // set minima and maxima (can be arbitrarily large)
    vector<double> tbqs_minima = { 0.0,          -TBQS_INFINITY, -TBQS_INFINITY, -TBQS_INFINITY };
    vector<double> tbqs_maxima = { TBQS_INFINITY, TBQS_INFINITY,  TBQS_INFINITY,  TBQS_INFINITY };

    // add EoS to vector
    chosen_EOSs.push_back( std::make_shared<EoS_conformal_cosh>(
                            c, T0, muB0, muS0, muQ0,
                            tbqs_minima, tbqs_maxima ) );
    init_conformal_diagonal(c, T0, muB0, muQ0, muS0);
    default_eos_name = "cosh_conformal";
  } else if ( settingsPtr->eos_type == "conformal_diagonal" )
  {
    const double Nc = 3.0, Nf = 2.5;  // u+d massless, s 'half massless'
    double c  = pi*pi*(2.0*(Nc*Nc-1.0)+(7.0/2.0)*Nc*Nf)/90.0;
    double T0 = 1.0, muB0 = 1.0, muQ0 = 1.0, muS0 = 1.0; // trivial scales
    init_conformal_diagonal(c, T0, muB0, muQ0, muS0);
    default_eos_name = "conformal_diagonal";
  } else if (settingsPtr->eos_type == "tanh_conformal")
  {
    const double Nc = 3.0, Nf = 2.5;  // u+d massless, s 'half massless'
    double c  = pi*pi*(2.0*(Nc*Nc-1.0)+(7.0/2.0)*Nc*Nf)/90.0;
    double T0 = 1.0, muB0 = 1.0, muQ0 = 1.0, muS0 = 1.0; // trivial scales
    double Tc = 220.0/constants::hbarc_MeVfm;
    double Ts = 120.0/constants::hbarc_MeVfm;
    init_tanh_conformal(c, T0, muB0, muQ0, muS0, Tc, Ts);
    default_eos_name = "tanh_conformal";
  } else if ( settingsPtr->eos_type == "table" )
  {
    eos_path = settingsPtr->eos_path;
    // add EoS to vector
    chosen_EOSs.push_back( std::make_shared<EoS_table>( eos_path ) );
    // setup tanh as first fallback for table EoS
    default_eos_name = "table";
    if ( use_tanh_conformal )
    {
      // pointer to default EoS (first element added above)
      pEoS_base p_default_EoS = chosen_EOSs.front();

      // look up grid maxima (without any extensions)
      std::vector<double> maxima = p_default_EoS->get_tbqs_maxima_no_ext();
      double Tmax   = maxima[0];
      double muBmax = maxima[1];
      double muQmax = maxima[2];
      double muSmax = maxima[3];

      // set overall scale using (Tmax,0,0,0)
      tbqs( Tmax, 0.0, 0.0, 0.0, p_default_EoS );
      double pTmax = pVal;
      double c     = pTmax / (Tmax*Tmax*Tmax*Tmax);

      // T-scale T0 = 1 by definition
      double T0 = 1.0;

      // set muB scale using (Tmax,muBmax,0,0)
      tbqs( Tmax, muBmax, 0.0, 0.0, p_default_EoS );
      //cout << pTmax << "   " << pVal << "   " << c << "   " << muBmax << endl;
      double muB0 = pow(c,0.25) * muBmax / sqrt( sqrt(pVal) - sqrt(pTmax) );

      // set muQ scale using (Tmax,0,muQmax,0)
      tbqs( Tmax, 0.0, muQmax, 0.0, p_default_EoS );
      //cout << pTmax << "   " << pVal << "   " << c << "   " << muQmax << endl;
      double muQ0 = pow(c,0.25) * muQmax / sqrt( sqrt(pVal) - sqrt(pTmax) );

      // set muS scale using (Tmax,0,0,muSmax)
      tbqs( Tmax, 0.0, 0.0, muSmax, p_default_EoS );
      //cout << pTmax << "   " << pVal << "   " << c << "   " << muSmax << endl;
      double muS0 = pow(c,0.25) * muSmax / sqrt( sqrt(pVal) - sqrt(pTmax) );

      // try rough scales for now (estimated by eye, not rigorously)
      double Tc = 220.0/constants::hbarc_MeVfm;
      double Ts = 120.0/constants::hbarc_MeVfm;

      init_tanh_conformal(c, T0, muB0, muQ0, muS0, Tc, Ts);
    }
    // use conformal as next fallback
    {
      formatted_output::update("Setting conformal equation of state as fallback");
      formatted_output::detail("all coefficients matched to p/T^4 at grid limits");

      // pointer to default EoS (first element added above)
      pEoS_base p_default_EoS = chosen_EOSs.front();

      // look up grid maxima (without any extensions)
      std::vector<double> maxima = p_default_EoS->get_tbqs_maxima_no_ext();
      double muBmax =    maxima[1];
      double muQmax =    maxima[2];
      double muSmax =    maxima[3]; 

	    double TmaxIni   = settingsPtr->Freeze_Out_Temperature/constants::hbarc_MeVfm;
	    double Tmax      = 1.1 * TmaxIni;
      //Tmax   = maxima[0];
      formatted_output::detail("Tmax = " + to_string(Tmax));
      // set overall scale using (Tmax,0,0,0)
      tbqs( Tmax, 0.0, 0.0, 0.0, p_default_EoS );
      double pTmax = pVal;
      double c     = pTmax / (Tmax*Tmax*Tmax*Tmax);

      // T-scale T0 = 1 by definition
      double T0 = 1.0;

      // set muB scale using (Tmax,muBmax,0,0)
      tbqs( Tmax, muBmax, 0.0, 0.0, p_default_EoS );
      //cout << pTmax << "   " << pVal << "   " << c << "   " << muBmax << endl;
      double muB0 = pow(c,0.25) * muBmax / sqrt( sqrt(pVal) - sqrt(pTmax) );

      // set muQ scale using (Tmax,0,muQmax,0)
      tbqs( Tmax, 0.0, muQmax, 0.0, p_default_EoS );
      //cout << pTmax << "   " << pVal << "   " << c << "   " << muQmax << endl;
      double muQ0 = pow(c,0.25) * muQmax / sqrt( sqrt(pVal) - sqrt(pTmax) );

      // set muS scale using (Tmax,0,0,muSmax)
      tbqs( Tmax, 0.0, 0.0, muSmax, p_default_EoS );
      //cout << pTmax << "   " << pVal << "   " << c << "   " << muSmax << endl;
      double muS0 = pow(c,0.25) * muSmax / sqrt( sqrt(pVal) - sqrt(pTmax) );
     init_conformal(c, T0, muB0, muQ0, muS0);
    }
    { // use diagonal conformal as final fallback
      formatted_output::update("Setting diagonal conformal equation of state "
                               "as final fallback");
      formatted_output::detail("all coefficients matched to p/T^4 at grid limits");

      // pointer to default EoS (first element added above)
      pEoS_base p_default_EoS = chosen_EOSs.front();

      // look up grid maxima (without any extensions)
      std::vector<double> maxima = p_default_EoS->get_tbqs_maxima_no_ext();
      //double Tmax   = maxima[0];
      double muBmax =  maxima[1];
      double muQmax =  maxima[2];
      double muSmax =  maxima[3]; 

	    double TmaxIni   = settingsPtr->Freeze_Out_Temperature/constants::hbarc_MeVfm;
	    double Tmax      = 1.1 * TmaxIni;
      //Tmax = maxima[0];

      // set overall scale using (Tmax,0,0,0)
      tbqs( Tmax, 0.0, 0.0, 0.0, p_default_EoS );
      double pTmax = pVal;
      double c  = pTmax / (Tmax*Tmax*Tmax*Tmax);

      // T-scale T0 = 1 by definition
      double T0 = 1.0;

      // set muB scale using (Tmax,muBmax,0,0)
      tbqs( Tmax, muBmax, 0.0, 0.0, p_default_EoS );
      //cout << pTmax << "   " << pVal << "   " << c << "   " << muBmax << endl;
      double muB0 = pow( c/(pVal - pTmax), 0.25) * muBmax;

      // set muQ scale using (Tmax,0,muQmax,0)
      tbqs( Tmax, 0.0, muQmax, 0.0, p_default_EoS );
      //cout << pTmax << "   " << pVal << "   " << c << "   " << muQmax << endl;
      double muQ0 = pow( c/(pVal - pTmax), 0.25) * muQmax;

      // set muS scale using (Tmax,0,0,muSmax)
      tbqs( Tmax, 0.0, 0.0, muSmax, p_default_EoS );
      //cout << pTmax << "   " << pVal << "   " << c << "   " << muSmax << endl;
      double muS0 = pow( c/(pVal - pTmax), 0.25) * muSmax;

      init_conformal_diagonal(c, T0, muB0, muQ0, muS0);

      default_eos_name = "table";

    }
  } else { // Unsupported EOS type - throw exception
    throw std::runtime_error("Unsupported equation of state type: " + settingsPtr->eos_type);
  }

  //============================================================================
  // create a map to access all needed EoSs by name
  // (this step *MUST BE DONE AFTER* chosen EoSs have been set,
  //  and each EoS must have a *UNIQUE NAME*)
  //============================================================================
  formatted_output::update( "Check order of equations of state:" );
  for ( auto & chosen_eos : chosen_EOSs )
  {
    formatted_output::detail( chosen_eos->name );
    chosen_EOS_map.insert({{ chosen_eos->name, chosen_eos }});
  }

	return;
}

 void EquationOfState::run_closure_test()
{
  const double hc = constants::hbarc_MeVfm;

  //==========================================================================
 /* std::cout << "Check conformal EoS:" << std::endl;
  for (double T0   = 0.0;     T0   <= 1200.01; T0   += 1200.0)
  for (double muB0 = -450.0; muB0 <= 450.01; muB0 += 450.0)
  for (double muS0 = -450.0; muS0 <= 450.01; muS0 += 450.0)
  for (double muQ0 = -450.0; muQ0 <= 450.01; muQ0 += 450.0)
  {
    std::vector<double> point = {T0/hc, muB0/hc, muQ0/hc, muS0/hc};
    std::vector<double> v = get_thermodynamics( point, "conformal" );
    std::cout << "Check conformal: " << T0 << "   " << muB0 << "   "
              << muQ0 << "   "<< muS0 << "   "
              << v[0]*hc*hc*hc*hc/(T0*T0*T0*T0) << "   "
              << v[1]*hc*hc*hc/(T0*T0*T0) << "   "
              << v[2]*hc*hc*hc/(T0*T0*T0) << "   "
              << v[3]*hc*hc*hc/(T0*T0*T0) << "   "
              << v[4]*hc*hc*hc/(T0*T0*T0) << "   "
              << v[5]*hc*hc*hc*hc/(T0*T0*T0*T0) << "   "
              << v[6] << "   "
              << v[7]*hc*hc/(T0*T0) << "   "
              << v[8]*hc*hc/(T0*T0) << "   "
              << v[9]*hc*hc/(T0*T0) << "   "
              << v[10]*hc*hc/(T0*T0) << "   "
              << v[11]*hc*hc/(T0*T0) << "   "
              << v[12]*hc*hc/(T0*T0) << "   "
              << v[13]*hc*hc/(T0*T0) << "   "
              << v[14]*hc*hc/(T0*T0) << "   "
              << v[15]*hc*hc/(T0*T0) << "   "
              << v[16]*hc*hc/(T0*T0) << std::endl;
  }
  std::cout << std::endl << std::endl << std::endl; */



/*
  //==========================================================================
  std::cout << "Check conformal_diagonal EoS:" << std::endl;
  for (double T0   = 0.0;     T0   <= 1200.01; T0   += 1200.0)
  for (double muB0 = -450.0; muB0 <= 450.01; muB0 += 450.0)
  for (double muS0 = -450.0; muS0 <= 450.01; muS0 += 450.0)
  for (double muQ0 = -450.0; muQ0 <= 450.01; muQ0 += 450.0)
  {
    std::vector<double> point = {T0/hc, muB0/hc, muQ0/hc, muS0/hc};
    std::vector<double> v = get_thermodynamics( point, "conformal_diagonal" );
    std::cout << "Check conformal_diagonal: "
              << T0 << "   " << muB0 << "   "
              << muQ0 << "   " << muS0 << "   "
              << v[0]*hc*hc*hc*hc/(T0*T0*T0*T0) << "   "
              << v[1]*hc*hc*hc/(T0*T0*T0) << "   "
              << v[2]*hc*hc*hc/(T0*T0*T0) << "   "
              << v[3]*hc*hc*hc/(T0*T0*T0) << "   "
              << v[4]*hc*hc*hc/(T0*T0*T0) << "   "
              << v[5]*hc*hc*hc*hc/(T0*T0*T0*T0) << "   "
              << v[6] << "   "
              << v[7]*hc*hc/(T0*T0) << "   "
              << v[8]*hc*hc/(T0*T0) << "   "
              << v[9]*hc*hc/(T0*T0) << "   "
              << v[10]*hc*hc/(T0*T0) << "   "
              << v[11]*hc*hc/(T0*T0) << "   "
              << v[12]*hc*hc/(T0*T0) << "   "
              << v[13]*hc*hc/(T0*T0) << "   "
              << v[14]*hc*hc/(T0*T0) << "   "
              << v[15]*hc*hc/(T0*T0) << "   "
              << v[16]*hc*hc/(T0*T0) << std::endl;
  }
  std::cout << std::endl << std::endl << std::endl;
*/

/*
  //==========================================================================
  std::cout << "Check non-conformal extension of table EoS:" << std::endl;
  std::cout << "Made it here" << std::endl;
  double e_In, rhoB_In, rhoS_In, rhoQ_In;
  for (double T0 =  5000.0; T0 <= 5000.01; T0 += 500.0)
  for (double muB0 = 2000.0; muB0 <= 2000.01; muB0 += 500.0)
  for (double muQ0 = -3000.0; muQ0 <= -(3000.0-0.01); muQ0 += 500.0)
  for (double muS0 = -2000.0; muS0 <= -(2000.0-0.01); muS0 += 500.0)
  {
    std::cout << "GETTING THERMODYNAMICS" << std::endl;
    std::vector<double> point = {T0/hc, muB0/hc, muQ0/hc, muS0/hc};
    std::vector<double> v = get_thermodynamics( point, "table" );
    std::cout << "Check nc_ext_table: " << T0 << "   " << muB0 << "   "
              << muQ0 << "   "<< muS0 << "   " << v[0] << "   " << v[6] << "   "
              << v[0]*hc*hc*hc*hc/(T0*T0*T0*T0) << std::endl;
    e_In    = v[5];
    rhoB_In = v[2];
    rhoS_In = v[3];
    rhoQ_In = v[4];
    std::cout << "GOT THERMODYNAMICS" << std::endl;
  }
  std::cout << std::endl << std::endl << std::endl;

  // closure test
  bool solution_found = false;
  double sLocal = s_out( e_In, rhoB_In, rhoS_In, rhoQ_In, solution_found );
  if ( solution_found )
    cout << "Closure test: successful!" << endl;
  else
    cout << "Closure test: unsuccessful!" << endl;

  cout << sLocal << "   " << e_In*hc << "   " << rhoB_In << "   "
        << rhoS_In << "   " << rhoQ_In << endl;

cout << "THERMO DUMP: "
    << tbqsPosition[0]*hc << "   " << tbqsPosition[1]*hc << "   "
    << tbqsPosition[2]*hc << "   " << tbqsPosition[3]*hc << "   "
    << pVal*hc << "   " << entrVal << "   " << BVal << "   "
    << SVal << "   " << QVal << "   " << eVal*hc << "   " << cs2Val << "   "
    << db2 << "   " << dq2 << "   " << ds2 << "   " << dbdq << "   "
    << dbds << "   " << dsdq << "   " << dtdb << "   " << dtdq << "   "
    << dtds << "   " << dt2 << endl;
*/

cout << "================================================================================\n"
      << "================================================================================\n"
      << "================================================================================\n";

  double Tstart   = 142.707028;
  double muBstart = -432.822491;
  double muSstart = -131.874802;
  double muQstart = 127.665770;
  double Tend     = Tstart;
  double muBend   = muBstart;
  double muSend   = muSstart;
  double muQend   = muQstart;

  for (double T0   = Tstart;   T0   <= Tend   + 0.01; T0   += 500.0)
  for (double muB0 = muBstart; muB0 <= muBend + 0.01; muB0 += 500.0)
  for (double muQ0 = muQstart; muQ0 <= muQend + 0.01; muQ0 += 500.0)
  for (double muS0 = muSstart; muS0 <= muSend + 0.01; muS0 += 500.0)
  {
    std::cout << "GETTING THERMODYNAMICS" << std::endl;
    std::vector<double> point = {T0/hc, muB0/hc, muQ0/hc, muS0/hc};

    // shared_ptr to EoS_table object
    pEoS_base table_EoS_object = chosen_EOS_map["table"];

    // call with debugging on (uses static library)
    std::dynamic_pointer_cast<EoS_table>(table_EoS_object)->set_debug_mode(true);
    std::vector<double> v = get_thermodynamics( point, "table" );
//    std::cout << "Check exact: " << T0 << "   " << muB0 << "   "
//              << muQ0 << "   "<< muS0 << "   " << v[5] << "   "
//              << v[2] << "   " << v[3] << "   " << v[4] << "   " << v[6] << std::endl;
    std::cout << "Check exact:";
    for (auto&e:v) cout << " " << e;
    std::cout << std::endl;

    // call with debugging on (uses interpolator)
    std::dynamic_pointer_cast<EoS_table>(table_EoS_object)->set_debug_mode(false);
    v = get_thermodynamics( point, "table" );
//    std::cout << "Check interpolant: " << T0 << "   " << muB0 << "   "
//              << muQ0 << "   "<< muS0 << "   " << v[5] << "   "
//              << v[2] << "   " << v[3] << "   " << v[4] << "   " << v[6] << std::endl;
    std::cout << "Check interpolant:";
    for (auto&e:v) cout << " " << e;
    std::cout << std::endl;
    
    std::cout << "GOT THERMODYNAMICS" << std::endl;
  }
  std::cout << std::endl << std::endl << std::endl;

}