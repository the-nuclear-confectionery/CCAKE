#include "eos.h"

//#include "read_in_hdf/read_in_hdf.h"
//#include "Stopwatch.h"
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

// functions calls to static EoS C library
//#include <lib.h>
#include "eos_delaunay/eos_delaunay.h"

#include "constants.h"
#include "eos_conformal.h"
#include "eos_extension.h"
#include "rootfinder.h"

using namespace constants;

using std::vector;
using std::string;

void EquationOfState::set_SettingsPtr( Settings * settingsPtr_in ) { settingsPtr = settingsPtr_in; }


void EquationOfState::init()
{
  cout << "Attempting read in of EoS from "
        << quantity_file << " and " << deriv_file << endl;
  init( quantity_file, deriv_file );
}

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
    vector<double> tbqs_minima = { 0.0,     -INFINITY, -INFINITY, -INFINITY };
    vector<double> tbqs_maxima = { INFINITY, INFINITY,  INFINITY,  INFINITY };

    // add EoS to vector
    chosen_EOSs.push_back( std::make_shared<EoS_conformal>(
                            c, T0, muB0, muS0, muQ0,
                            tbqs_minima, tbqs_maxima ) );
  }
  else if ( settingsPtr->EoS_type == "Table" )
  {
    // add EoS to vector
    chosen_EOSs.push_back( std::make_shared<EoS_table>( quantityFile, derivFile ) );
  }
  //////////////////////////////////////////////////////////////////////////////



  //////////////////////////////////////////////////////////////////////////////
  // SET REMAINING BACKUP EQUATIONS OF STATE
  // - non-conformal extension (extended everywhere)
  // - conformal extension (extended everywhere)
  // - purely conformal fallback (always use this to guarantee solution)
  //////////////////////////////////////////////////////////////////////////////
  if ( settingsPtr->EoS_type != "Conformal" ) // use conformal as fallback
  {
    std::cout << "Setting conformal equation of state as fallback" << std::endl;
    const double Nc = 3.0, Nf = 2.5;  // u+d massless, s 'half massless'
    double c  = pi*pi*(2.0*(Nc*Nc-1.0)+(7.0/2.0)*Nc*Nf)/90.0;
    double T0 = 1.0, muB0 = 1.0, muQ0 = 1.0, muS0 = 1.0; // trivial scales

    // set minima and maxima (can be arbitrarily large)
    vector<double> tbqs_minima = { 0.0,     -INFINITY, -INFINITY, -INFINITY };
    vector<double> tbqs_maxima = { INFINITY, INFINITY,  INFINITY,  INFINITY };

    // add EoS to vector
    chosen_EOSs.push_back( std::make_shared<EoS_conformal>(
                            c, T0, muB0, muS0, muQ0,
                            tbqs_minima, tbqs_maxima ) );
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

