#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "constants.h"
#include "eos.h"
#include "eos_base.h"
#include "eos_conformal.h"
#include "eos_header.h"
#include "eos_nonconformal_extension.h"
#include "eos_table.h"
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
  {
    const double hc = constants::hbarc_MeVfm;

    //==========================================================================
    std::cout << "Check conformal EoS:" << std::endl;
    for (double T0   = 0.0;     T0   <= 1000.01; T0   += 1000.0)
    for (double muB0 = -1000.0; muB0 <= 1000.01; muB0 += 1000.0)
    for (double muS0 = -1000.0; muS0 <= 1000.01; muS0 += 1000.0)
    for (double muQ0 = -1000.0; muQ0 <= 1000.01; muQ0 += 1000.0)
    {
      std::vector<double> point = {T0/hc, muB0/hc, muQ0/hc, muS0/hc};
      std::vector<double> v = get_thermodynamics( point, "conformal" );
      std::cout << "Check conformal: " << T0 << "   " << muB0 << "   "
                << muQ0 << "   "<< muS0 << "   " << v[0] << std::endl;
    }
    std::cout << std::endl << std::endl << std::endl;

    //==========================================================================
    std::cout << "Check non-conformal extension of table EoS:" << std::endl;
    for (double T0 =  75.0; T0 <= 75.01; T0 += 500.0)
    for (double muB0 = 1000.0; muB0 <= 1000.01; muB0 += 500.0)
    for (double muS0 = 2000.0; muS0 <= 2000.01; muS0 += 500.0)
    for (double muQ0 = 3000.0; muQ0 <= 3000.01; muQ0 += 500.0)
//    for (double T0 =  150.0-0.0000001; T0 <= 150.01; T0 += 500.0)
//    for (double muB0 = 150.0-0.0000001; muB0 <= 150.01; muB0 += 500.0)
//    for (double muS0 = 300.0-0.0000001; muS0 <= 300.01; muS0 += 500.0)
//    for (double muQ0 = 450.0-0.0000001; muQ0 <= 450.01; muQ0 += 500.0)
//    for (double T0 =  0.0; T0 <= 100000.01; T0 += 500.0)
//    for (double muB0 = 0.0; muB0 <= 0.01; muB0 += 500.0)
//    for (double muS0 = 0.0; muS0 <= 0.01; muS0 += 500.0)
//    for (double muQ0 = 0.0; muQ0 <= 0.01; muQ0 += 500.0)
    {
      std::vector<double> point = {T0/hc, muB0/hc, muQ0/hc, muS0/hc};
      std::vector<double> v = get_thermodynamics( point, "table" );
      std::cout << "Check nc_ext_table: " << T0 << "   " << muB0 << "   "
                << muQ0 << "   "<< muS0 << "   " << v[0] << "   " << v[6] << "   "
                << v[0]*hc*hc*hc*hc/(T0*T0*T0*T0) << std::endl;
    }
    std::cout << std::endl << std::endl << std::endl;

    exit(11);
  }
}

////////////////////////////////////////////////////////////////////////////////
void EquationOfState::init(string quantityFile, string derivFile)
{
	tbqsPosition.resize(4);

  //////////////////////////////////////////////////////////////////////////////
  // THIS IF-CHAIN INITIALIZES THE DEFAULT EOS TO USE
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
  else if ( settingsPtr->EoS_type == "Houston" )
  {
    // add EoS to vector
    chosen_EOSs.push_back( std::make_shared<EoS_table>( quantityFile, derivFile ) );
    default_eos_name = "table";
  }
  //////////////////////////////////////////////////////////////////////////////



  //////////////////////////////////////////////////////////////////////////////
  // SET REMAINING BACKUP EQUATIONS OF STATE
  // - non-conformal extension (extended everywhere)
  // - conformal extension (extended everywhere)
  // - purely conformal fallback (always use this to guarantee solution)
  //////////////////////////////////////////////////////////////////////////////

  // set non-conformal extension
  if ( settingsPtr->EoS_type != "Conformal" ) // redundant with below, but oh well
  {
    // pointer to default EoS (first element added above)
    pEoS_base p_default_EoS = chosen_EOSs.front();
    
    // set minima and maxima (can be arbitrarily large)
    vector<double> tbqs_minima = { 0.0,          -TBQS_INFINITY, -TBQS_INFINITY, -TBQS_INFINITY };
    vector<double> tbqs_maxima = { TBQS_INFINITY, TBQS_INFINITY,  TBQS_INFINITY,  TBQS_INFINITY };

    chosen_EOSs.push_back(
        std::make_shared<EoS_nonconformal_extension>(
                          p_default_EoS, tbqs_minima, tbqs_maxima ) );
  }


  // use conformal as fallback
  if ( settingsPtr->EoS_type != "Conformal" )
  {
    std::cout << "Setting conformal equation of state as fallback" << std::endl;
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
  }




  //////////////////////////////////////////////////////////////////////////////
  // create a map to access all needed EoSs by name
  // (this step *MUST BE DONE AFTER* chosen EoSs have been set,
  //  and each EoS must have a *UNIQUE NAME*)
  for ( auto & chosen_eos : chosen_EOSs )
  {
    std::cout << "Before " << chosen_eos->name << ": chosen_EOS_map.size = " << chosen_EOS_map.size() << std::endl;
    chosen_EOS_map.insert({{ chosen_eos->name, chosen_eos }});
    std::cout << "After " << chosen_eos->name << ": chosen_EOS_map.size = " << chosen_EOS_map.size() << std::endl;
  }

  




  //////////////////////////////////////////////////////////////////////////////
  // set method for locating point in phase diagram
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
  //////////////////////////////////////////////////////////////////////////////

	return;
}

