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
#include <lib.h>
#include "eos_delaunay/eos_delaunay.h"

#include "constants.h"
#include "eos_conformal.h"
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
  if ( settingsPtr->EoS_type == "Conformal" )
  {
    std::cout << "Setting up equation of state for Gubser checks" << std::endl;
    const double Nc = 3.0, Nf = 2.5;  // u+d massless, s 'half massless'
    eos_conformal::c    = pi*pi*(2.0*(Nc*Nc-1.0)+(7.0/2.0)*Nc*Nf)/90.0;
    eos_conformal::T0   = 1.0;
    eos_conformal::muB0 = 1.0;
    eos_conformal::muQ0 = 1.0;
    eos_conformal::muS0 = 1.0;
    
    // set density-computing functions appropriately
    std::function<void(double[], double[])> f_eBSQ = eos_conformal::get_eBSQ;
    set_conformal_eBSQ_functional( f_eBSQ );
    std::function<void(double[], double[])> f_sBSQ = eos_conformal::get_sBSQ;
    set_conformal_sBSQ_functional( f_sBSQ );

    /*
    // set minima and maxima (can be arbitrarily large)
    minT   =  0.0;     minMuB = -10000.0;
    minMuS = -10000.0; minMuQ = -10000.0;
    maxT   =  10000.0; maxMuB =  10000.0;
    maxMuS =  10000.0; maxMuQ =  10000.0;
    */
    double conformal_scale = 100000.0;
    conformal_tbqs_minima = { 0.0,             -conformal_scale,
                             -conformal_scale, -conformal_scale };
    conformal_tbqs_maxima = { conformal_scale, conformal_scale,
                              conformal_scale, conformal_scale };
    
  }
  else if ( use_static_C_library )
  {
    // initialize things needed to use static C library
    cout << "Initializing EoS C library" << endl;
    //initialize("/projects/jnorhos/BSQ/EoS_BQS_Derivatives/Coefficients_Parameters.dat");
    initialize("../EoS_BQS_Derivatives/Coefficients_Parameters.dat");

    // set density-computing functions appropriately
    std::function<void(double[], double[])> f_eBSQ = STANDARD_get_eBSQ_densities;
    set_eBSQ_functional( f_eBSQ );
    std::function<void(double[], double[])> f_sBSQ = STANDARD_get_sBSQ_densities;
    set_sBSQ_functional( f_sBSQ );

    // load EoS tables, assess grid range
    std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    init_grid_ranges_only(quantityFile, derivFile);
  }
  else  // if thermo not from static C library, read in table from file
  {
    // initialize things needed to use static C library
    equation_of_state_table_filename = "./EoS/Houston/Default/thermo.dat";
    cout << "Initializing EoS from input file(s): "
        << equation_of_state_table_filename << endl;
    equation_of_state_table.initialize( equation_of_state_table_filename );

    // set names of EoS quantities to interpolate, in order
    equation_of_state_table.set_grid_names(
      vector<string>{ "T","muB","muQ","muS" } );
    equation_of_state_table.set_field_names(
      vector<string>{ "p","s","B","S","Q","e","cs2",
                      "chiBB","chiQQ","chiSS","chiBQ","chiBS",
                      "chiQS","chiTB","chiTQ","chiTS","chiTT" } );

    // finally, all dimensionalful quantities should be convert to fm and
    // all dimensionless quantities need to be rescaled appropriately
    /*equation_of_state_table.rescale_axis( "T",   1.0/hbarc_MeVfm );
    equation_of_state_table.rescale_axis( "muB", 1.0/hbarc_MeVfm );
    equation_of_state_table.rescale_axis( "muQ", 1.0/hbarc_MeVfm );
    equation_of_state_table.rescale_axis( "muS", 1.0/hbarc_MeVfm );*/
    equation_of_state_table.rescale_axes( 1.0/hbarc_MeVfm );  // do all axes at once

    equation_of_state_table.rescale( "p",     "T", 4 );
    equation_of_state_table.rescale( "e",     "T", 4 );
    equation_of_state_table.rescale( "s",     "T", 3 );
    equation_of_state_table.rescale( "B",     "T", 3 );
    equation_of_state_table.rescale( "S",     "T", 3 );
    equation_of_state_table.rescale( "Q",     "T", 3 );
    equation_of_state_table.rescale( "chiBB", "T", 2 );
    equation_of_state_table.rescale( "chiQQ", "T", 2 );
    equation_of_state_table.rescale( "chiSS", "T", 2 );
    equation_of_state_table.rescale( "chiBQ", "T", 2 );
    equation_of_state_table.rescale( "chiBS", "T", 2 );
    equation_of_state_table.rescale( "chiQS", "T", 2 );
    equation_of_state_table.rescale( "chiTB", "T", 2 );
    equation_of_state_table.rescale( "chiTQ", "T", 2 );
    equation_of_state_table.rescale( "chiTS", "T", 2 );
    equation_of_state_table.rescale( "chiTT", "T", 2 );

    // set grid ranges
    /*vector<double> grid_minima = equation_of_state_table.get_grid_minima();
    vector<double> grid_maxima = equation_of_state_table.get_grid_maxima();
    minT   = grid_minima[0]; minMuB = grid_minima[1];
    minMuS = grid_minima[2]; minMuQ = grid_minima[3];
    maxT   = grid_maxima[0]; maxMuB = grid_maxima[1];
    maxMuS = grid_maxima[2]; maxMuQ = grid_maxima[3];*/
    tbqs_minima = equation_of_state_table.get_grid_minima();
    tbqs_maxima = equation_of_state_table.get_grid_maxima();

    // set functions from interpolant
    std::function<void(double[], double[])> f_eBSQ = get_eBSQ_densities_from_interpolator;
    set_eBSQ_functional( f_eBSQ );
    std::function<void(double[], double[])> f_sBSQ = get_sBSQ_densities_from_interpolator;
    set_sBSQ_functional( f_sBSQ );

    // load EoS tables, assess grid range
    //std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    //init_grid(quantityFile, derivFile);
  }

  // set up conformal option in this case as well
  if ( use_conformal_as_fallback )
  {
    std::cout << "Setting up equation of state as fallback "
                 "when default equation of state fails" << std::endl;
    const double Nc = 3.0, Nf = 2.5;  // u+d massless, s 'half massless'
    eos_conformal::c    = pi*pi*(2.0*(Nc*Nc-1.0)+(7.0/2.0)*Nc*Nf)/90.0;
    eos_conformal::T0   = 1.0;
    eos_conformal::muB0 = 1.0;
    eos_conformal::muQ0 = 1.0;
    eos_conformal::muS0 = 1.0;
    
    // set density-computing functions appropriately
    std::function<void(double[], double[])> f_eBSQ = eos_conformal::get_eBSQ;
    set_conformal_eBSQ_functional( f_eBSQ );
    std::function<void(double[], double[])> f_sBSQ = eos_conformal::get_sBSQ;
    set_conformal_sBSQ_functional( f_sBSQ );

    // set minima and maxima (can be arbitrarily large)
    double conformal_scale = 100000.0;
    conformal_tbqs_minima = { 0.0,             -conformal_scale,
                             -conformal_scale, -conformal_scale };
    conformal_tbqs_maxima = { conformal_scale, conformal_scale,
                              conformal_scale, conformal_scale };
    
  }
  //////////////////////////////////////////////////////////////////////////////

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

void EquationOfState::init_grid_ranges_only(string quantityFile, string derivFile)
{
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    std::ifstream dataFile;
    dataFile.open(quantityFile);

    double maxMuB        = 0.0;
    double minMuB        = 0.0;
    double maxMuQ        = 0.0;
    double minMuQ        = 0.0;
    double maxMuS        = 0.0;
    double minMuS        = 0.0;
    double maxT          = 0.0;
    double minT          = 0.0;
    double tit, muBit, muQit, muSit, pit, entrit, bit, sit, qit, eit, cs2it;

    int count = 0;
    while ( dataFile >> tit >> muBit >> muQit >> muSit
                     >> pit >> entrit >> bit >> sit >> qit >> eit >> cs2it )
    {

		// Christopher Plumberg:
		// put T and mu_i in units of 1/fm
		tit   /= hbarc_MeVfm;
		muBit /= hbarc_MeVfm;
		muSit /= hbarc_MeVfm;
		muQit /= hbarc_MeVfm;

    if(count++ == 0)
    {
      minT   = tit;
      maxT   = tit;
      minMuB = muBit;
      maxMuB = muBit;     //initialize eos range variables
      minMuQ = muQit;
      maxMuQ = muQit;
      minMuS = muSit;
      maxMuS = muSit;
    }
        
		if (count%100000==0) std::cout << "Read in line# " << count << std::endl;
		
        if (maxT < tit) maxT = tit;
        if (minT > tit) minT = tit;
        if (maxMuB < muBit) maxMuB = muBit;
        if (minMuB > muBit) minMuB = muBit;
        if (maxMuQ < muQit) maxMuQ = muQit;
        if (minMuQ > muQit) minMuQ = muQit;
        if (maxMuS < muSit) maxMuS = muSit;
        if (minMuS > muSit) minMuS = muSit;
        
	}

  // Add this for now
  if (true)
  {
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << "!!           WARNING           !!" << endl;
    cout << "!!   FORCING EOS GRID RANGES   !!" << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    minMuB = -8000.0/hbarc_MeVfm;
    minMuQ = -8000.0/hbarc_MeVfm;
    minMuS = -8000.0/hbarc_MeVfm;
    maxMuB = 8000.0/hbarc_MeVfm;
    maxMuQ = 8000.0/hbarc_MeVfm;
    maxMuS = 8000.0/hbarc_MeVfm;
  }

  // initialize grid ranges here
  tbqs_minima = {minT, minMuB, minMuQ, minMuS};
  tbqs_maxima = {maxT, maxMuB, maxMuQ, maxMuS};

  dataFile.close();

	std::cout << "All initializations finished!" << std::endl;

  return;
}