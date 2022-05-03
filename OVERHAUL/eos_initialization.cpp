#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <iostream>
#include <memory>
#include <vector>
#include <fstream>
#include <string>

#include "constants.h"
#include "eos.h"
#include "eos_base.h"
#include "eos_conformal.h"
#include "eos_conformal_diagonal.h"
#include "eos_header.h"
#include "eos_nonconformal_extension.h"
#include "eos_table.h"
#include "eos_tanh_conformal.h"
#include "eos_delaunay/eos_delaunay.h"
#include "rootfinder.h"

using namespace constants;

using std::vector;
using std::string;

////////////////////////////////////////////////////////////////////////////////
void EquationOfState::set_SettingsPtr( Settings * settingsPtr_in ) { settingsPtr = settingsPtr_in; }

////////////////////////////////////////////////////////////////////////////////
void EquationOfState::init()
{
  cout << "Attempting read in of EoS from "
        << quantity_file << " and " << deriv_file << endl;
  init( quantity_file, deriv_file );

  bool do_eos_checks = true;
  if ( do_eos_checks )
    run_closure_test();
}

////////////////////////////////////////////////////////////////////////////////
void EquationOfState::init(string quantityFile, string derivFile)
{
	tbqsPosition.resize(4);

  //============================================================================
  // THIS IF-CHAIN INITIALIZES THE DEFAULT EOS TO USE
  //============================================================================
  // SET UP CONFORMAL EOS
  if ( settingsPtr->EoS_type == "Conformal" )
  {
    std::cout << "Setting up equation of state for Gubser checks" << std::endl;
    const double Nc = 3.0, Nf = 2.5;  // u+d massless, s 'half massless'
    double c  = pi*pi*(2.0*(Nc*Nc-1.0)+(7.0/2.0)*Nc*Nf)/90.0;
    double T0 = 1.0, muB0 = 1.0, muQ0 = 1.0, muS0 = 1.0; // trivial scales

    // set minima and maxima (can be arbitrarily large)
    vector<double> tbqs_minima = { 0.0,          -TBQS_INFINITY, -TBQS_INFINITY, -TBQS_INFINITY };
    vector<double> tbqs_maxima = { TBQS_INFINITY, TBQS_INFINITY,  TBQS_INFINITY,  TBQS_INFINITY };

    // add EoS to vector
    chosen_EOSs.push_back( std::make_shared<EoS_conformal>(
                            c, T0, muB0, muS0, muQ0,
                            tbqs_minima, tbqs_maxima ) );
    default_eos_name = "conformal";
  }
  // SET UP CONFORMAL DIAGONAL EOS
  if ( settingsPtr->EoS_type == "Conformal_Diagonal" )
  {
    std::cout << "Setting DIAGONAL conformal equation of state" << std::endl;
    const double Nc = 3.0, Nf = 2.5;  // u+d massless, s 'half massless'
    double c  = pi*pi*(2.0*(Nc*Nc-1.0)+(7.0/2.0)*Nc*Nf)/90.0;
    double T0 = 1.0, muB0 = 1.0, muQ0 = 1.0, muS0 = 1.0; // trivial scales

    // set minima and maxima for rootfinder (can be arbitrarily large)
    vector<double> tbqs_minima = { 0.0,          -TBQS_INFINITY, -TBQS_INFINITY, -TBQS_INFINITY };
    vector<double> tbqs_maxima = { TBQS_INFINITY, TBQS_INFINITY,  TBQS_INFINITY,  TBQS_INFINITY };

    cout << "DIAGONAL Conformal fallback EoS set up with following parameters:" << endl;
    cout << "  --> c    = " << c << endl;
    cout << "  --> T0   = " << T0 << endl;
    cout << "  --> muB0 = " << muB0 << endl;
    cout << "  --> muQ0 = " << muQ0 << endl;
    cout << "  --> muS0 = " << muS0 << endl;

    // add matched conformal EoS to vector of EoSs
    chosen_EOSs.push_back( std::make_shared<EoS_conformal_diagonal>(
                            c, T0, muB0, muS0, muQ0,
                            tbqs_minima, tbqs_maxima, "conformal_diagonal" ) );

    conformal_diagonal_EoS = EoS_conformal_diagonal( c, T0, muB0, muS0, muQ0,
                                                     tbqs_minima, tbqs_maxima,
                                                     "conformal_diagonal" );
    default_eos_name = "conformal_diagonal";

  }
  // SET UP HOUSTON TABLE EOS
  else if ( settingsPtr->EoS_type == "Houston" )
  {
    // add EoS to vector
    chosen_EOSs.push_back( std::make_shared<EoS_table>( quantityFile, derivFile ) );
    default_eos_name = "table";
  }
  //============================================================================
  //============================================================================




  //============================================================================
  // FOR TABLE/HOUSTON EOS, define fallbacks if default fails
  //============================================================================
  if ( settingsPtr->EoS_type == "Houston" || settingsPtr->EoS_type == "Table" )
  {
    //==========================================================================
    // use tanh-modulated "conformal" as first fallback for table/Houston EoS
    //==========================================================================
    if ( use_tanh_conformal )
    {
      std::cout << "Setting tanh-modulated \"conformal\" equation of state as fallback" << std::endl;
      std::cout << "  --> all coefficients matched to p/T^4 at grid limits" << std::endl;

      // pointer to default EoS (first element added above)
      pEoS_base p_default_EoS = chosen_EOSs.front();

      // look up grid maxima (without any extensions)
      std::vector<double> maxima =  p_default_EoS->get_tbqs_maxima_no_ext();
      double Tmax = maxima[0];
      double muBmax = maxima[1];
      double muQmax = maxima[2];
      double muSmax = maxima[3];

      // set overall scale using (Tmax,0,0,0)
      tbqs( Tmax, 0.0, 0.0, 0.0, p_default_EoS );
      double pTmax = pVal;
      double c  = pTmax / (Tmax*Tmax*Tmax*Tmax);

      // T-scale T0 = 1 by definition
      double T0 = 1.0;

      // set muB scale using (Tmax,muBmax,0,0)
      tbqs( Tmax, muBmax, 0.0, 0.0, p_default_EoS );
      cout << pTmax << "   " << pVal << "   " << c << "   " << muBmax << endl;
      double muB0 = pow(c,0.25) * muBmax / sqrt( sqrt(pVal) - sqrt(pTmax) );

      // set muQ scale using (Tmax,0,muQmax,0)
      tbqs( Tmax, 0.0, muQmax, 0.0, p_default_EoS );
      cout << pTmax << "   " << pVal << "   " << c << "   " << muQmax << endl;
      double muQ0 = pow(c,0.25) * muQmax / sqrt( sqrt(pVal) - sqrt(pTmax) );

      // set muS scale using (Tmax,0,0,muSmax)
      tbqs( Tmax, 0.0, 0.0, muSmax, p_default_EoS );
      cout << pTmax << "   " << pVal << "   " << c << "   " << muSmax << endl;
      double muS0 = pow(c,0.25) * muSmax / sqrt( sqrt(pVal) - sqrt(pTmax) );

      // try rough scales for now (estimated by eye, not rigorously)
      double Tc = 220.0/constants::hbarc_MeVfm;
      double Ts = 120.0/constants::hbarc_MeVfm;

      // set minima and maxima for rootfinder (can be arbitrarily large)
      vector<double> tbqs_minima = { 0.0,          -TBQS_INFINITY, -TBQS_INFINITY, -TBQS_INFINITY };
      vector<double> tbqs_maxima = { TBQS_INFINITY, TBQS_INFINITY,  TBQS_INFINITY,  TBQS_INFINITY };

      cout << "Tanh-modulated \"conformal\" fallback EoS set up with following parameters:" << endl;
      cout << "  --> c    = " << c << endl;
      cout << "  --> T0   = " << T0 << endl;
      cout << "  --> muB0 = " << muB0 << endl;
      cout << "  --> muQ0 = " << muQ0 << endl;
      cout << "  --> muS0 = " << muS0 << endl;
      cout << "  --> Tc   = " << Tc << endl;
      cout << "  --> Ts   = " << Ts << endl;

      // add matched conformal EoS to vector of EoSs
      chosen_EOSs.push_back( std::make_shared<EoS_tanh_conformal>(
                              c, T0, muB0, muS0, muQ0, Tc, Ts,
                              tbqs_minima, tbqs_maxima, "tanh_conformal" ) );
    }



    //==========================================================================
    // use conformal as next fallback
    //==========================================================================
    {
      std::cout << "Setting conformal equation of state as fallback" << std::endl;
      std::cout << "  --> all coefficients matched to p/T^4 at grid limits" << std::endl;

      // pointer to default EoS (first element added above)
      pEoS_base p_default_EoS = chosen_EOSs.front();

      // look up grid maxima (without any extensions)
      std::vector<double> maxima =  p_default_EoS->get_tbqs_maxima_no_ext();
      double Tmax = maxima[0];
      double muBmax = maxima[1];
      double muQmax = maxima[2];
      double muSmax = maxima[3];

      // set overall scale using (Tmax,0,0,0)
      tbqs( Tmax, 0.0, 0.0, 0.0, p_default_EoS );
      double pTmax = pVal;
      double c  = pTmax / (Tmax*Tmax*Tmax*Tmax);

      //const double hc = constants::hbarc_MeVfm;

      // T-scale T0 = 1 by definition
      double T0 = 1.0;

      // set muB scale using (Tmax,muBmax,0,0)
      tbqs( Tmax, muBmax, 0.0, 0.0, p_default_EoS );
      cout << pTmax << "   " << pVal << "   " << c << "   " << muBmax << endl;
      double muB0 = pow(c,0.25) * muBmax / sqrt( sqrt(pVal) - sqrt(pTmax) );

      // set muQ scale using (Tmax,0,muQmax,0)
      tbqs( Tmax, 0.0, muQmax, 0.0, p_default_EoS );
      cout << pTmax << "   " << pVal << "   " << c << "   " << muQmax << endl;
      double muQ0 = pow(c,0.25) * muQmax / sqrt( sqrt(pVal) - sqrt(pTmax) );

      // set muS scale using (Tmax,0,0,muSmax)
      tbqs( Tmax, 0.0, 0.0, muSmax, p_default_EoS );
      cout << pTmax << "   " << pVal << "   " << c << "   " << muSmax << endl;
      double muS0 = pow(c,0.25) * muSmax / sqrt( sqrt(pVal) - sqrt(pTmax) );

      // set minima and maxima for rootfinder (can be arbitrarily large)
      vector<double> tbqs_minima = { 0.0,          -TBQS_INFINITY, -TBQS_INFINITY, -TBQS_INFINITY };
      vector<double> tbqs_maxima = { TBQS_INFINITY, TBQS_INFINITY,  TBQS_INFINITY,  TBQS_INFINITY };

      cout << "Conformal fallback EoS set up with following parameters:" << endl;
      cout << "  --> c    = " << c << endl;
      cout << "  --> T0   = " << T0 << endl;
      cout << "  --> muB0 = " << muB0 << endl;
      cout << "  --> muQ0 = " << muQ0 << endl;
      cout << "  --> muS0 = " << muS0 << endl;

      // add matched conformal EoS to vector of EoSs
      chosen_EOSs.push_back( std::make_shared<EoS_conformal>(
                              c, T0, muB0, muS0, muQ0,
                              tbqs_minima, tbqs_maxima, "conformal" ) );
    }

  }




    //==========================================================================
    // use diagonal conformal as final fallback (MUST ALWAYS INCLUDE)
    //==========================================================================
    if ( settingsPtr->EoS_type != "Conformal_Diagonal" )
    {
      std::cout << "Setting DIAGONAL conformal equation of state as FINAL fallback" << std::endl;
      std::cout << "  --> all coefficients matched to p/T^4 at grid limits" << std::endl;

      // pointer to default EoS (first element added above)
      pEoS_base p_default_EoS = chosen_EOSs.front();

      // look up grid maxima (without any extensions)
      std::vector<double> maxima =  p_default_EoS->get_tbqs_maxima_no_ext();
      double Tmax = maxima[0];
      double muBmax = maxima[1];
      double muQmax = maxima[2];
      double muSmax = maxima[3];

      // set overall scale using (Tmax,0,0,0)
      tbqs( Tmax, 0.0, 0.0, 0.0, p_default_EoS );
      double pTmax = pVal;
      double c  = pTmax / (Tmax*Tmax*Tmax*Tmax);

      //const double hc = constants::hbarc_MeVfm;

      // T-scale T0 = 1 by definition
      double T0 = 1.0;

      // set muB scale using (Tmax,muBmax,0,0)
      tbqs( Tmax, muBmax, 0.0, 0.0, p_default_EoS );
      cout << pTmax << "   " << pVal << "   " << c << "   " << muBmax << endl;
      double muB0 = pow( c/(pVal - pTmax), 0.25) * muBmax;

      // set muQ scale using (Tmax,0,muQmax,0)
      tbqs( Tmax, 0.0, muQmax, 0.0, p_default_EoS );
      cout << pTmax << "   " << pVal << "   " << c << "   " << muQmax << endl;
      double muQ0 = pow( c/(pVal - pTmax), 0.25) * muQmax;

      // set muS scale using (Tmax,0,0,muSmax)
      tbqs( Tmax, 0.0, 0.0, muSmax, p_default_EoS );
      cout << pTmax << "   " << pVal << "   " << c << "   " << muSmax << endl;
      double muS0 = pow( c/(pVal - pTmax), 0.25) * muSmax;

      // set minima and maxima for rootfinder (can be arbitrarily large)
      vector<double> tbqs_minima = { 0.0,          -TBQS_INFINITY, -TBQS_INFINITY, -TBQS_INFINITY };
      vector<double> tbqs_maxima = { TBQS_INFINITY, TBQS_INFINITY,  TBQS_INFINITY,  TBQS_INFINITY };

      cout << "Conformal fallback EoS set up with following parameters:" << endl;
      cout << "  --> c    = " << c << endl;
      cout << "  --> T0   = " << T0 << endl;
      cout << "  --> muB0 = " << muB0 << endl;
      cout << "  --> muQ0 = " << muQ0 << endl;
      cout << "  --> muS0 = " << muS0 << endl;

      // add matched conformal EoS to vector of EoSs
      chosen_EOSs.push_back( std::make_shared<EoS_conformal_diagonal>(
                              c, T0, muB0, muS0, muQ0,
                              tbqs_minima, tbqs_maxima, "conformal_diagonal" ) );

      conformal_diagonal_EoS = EoS_conformal_diagonal( c, T0, muB0, muS0, muQ0,
                                                       tbqs_minima, tbqs_maxima,
                                                       "conformal_diagonal" );
    }







  //============================================================================
  // create a map to access all needed EoSs by name
  // (this step *MUST BE DONE AFTER* chosen EoSs have been set,
  //  and each EoS must have a *UNIQUE NAME*)
  //============================================================================
  for ( auto & chosen_eos : chosen_EOSs )
  {
    std::cout << "Before " << chosen_eos->name
              << ": chosen_EOS_map.size = "
              << chosen_EOS_map.size() << std::endl;

    chosen_EOS_map.insert({{ chosen_eos->name, chosen_eos }});

    std::cout << "After " << chosen_eos->name
              << ": chosen_EOS_map.size = "
              << chosen_EOS_map.size() << std::endl;
  }

  




  //============================================================================
  // set method for locating point in phase diagram
  //============================================================================
  if ( use_rootfinder )
  {
    // initialize Rootfinder ranges
    //rootfinder.set_grid_ranges( minT, maxT, minMuB, maxMuB,
    //                            minMuS, maxMuS, minMuQ, maxMuQ );
  }
  else if ( use_delaunay )
  {
    // initialize corresponding interpolator for each table
    cout << "Initialize Delaunay interpolators" << endl;
    e_delaunay.init(    quantityFile, 0 );	// 0 - energy density
    entr_delaunay.init( quantityFile, 1 );	// 1 - entropy density
  }
  //============================================================================

	return;
}



void EquationOfState::run_closure_test()
{
  const double hc = constants::hbarc_MeVfm;
/*
  //==========================================================================
  std::cout << "Check conformal EoS:" << std::endl;
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
  std::cout << std::endl << std::endl << std::endl;
*/


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

  for (double T0 =  436.982; T0 <= 436.982 + 0.01; T0 += 500.0)
  for (double muB0 = 0.0; muB0 <= 0.01; muB0 += 500.0)
  for (double muQ0 = 0.0; muQ0 <= 0.01; muQ0 += 500.0)
  for (double muS0 = 0.0; muS0 <= 0.01; muS0 += 500.0)
  {
    std::cout << "GETTING THERMODYNAMICS" << std::endl;
    std::vector<double> point = {T0/hc, muB0/hc, muQ0/hc, muS0/hc};

    // shared_ptr to EoS_table object
    pEoS_base table_EoS_object = chosen_EOS_map["table"];

    // call with debugging on (uses static library)
    std::dynamic_pointer_cast<EoS_table>(table_EoS_object)->set_debug_mode(true);
    std::vector<double> v = get_thermodynamics( point, "table" );
    std::cout << "Check exact: " << T0 << "   " << muB0 << "   "
              << muQ0 << "   "<< muS0 << "   " << v[5] << "   "
              << v[2] << "   " << v[3] << "   " << v[4] << std::endl;

    // call with debugging on (uses interpolator)
    std::dynamic_pointer_cast<EoS_table>(table_EoS_object)->set_debug_mode(false);
    v = get_thermodynamics( point, "table" );
    std::cout << "Check interpolant: " << T0 << "   " << muB0 << "   "
              << muQ0 << "   "<< muS0 << "   " << v[5] << "   "
              << v[2] << "   " << v[3] << "   " << v[4] << std::endl;
    
    std::cout << "GOT THERMODYNAMICS" << std::endl;
  }
  std::cout << std::endl << std::endl << std::endl;


  exit(11);
}