#include "eos_table.h"

InterpolatorND<4> EoS_table::equation_of_state_table;

#include <algorithm>
#include <functional>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "eos_extension.h"
#include "eos_header.h"
#include "eos_table.h"

EoS_table::EoS_table( string quantityFile, string derivFile )
{
  // allow to use static C library instead of table
  if ( use_static_C_library )
  {
    // initialize things needed to use static C library
    cout << "Initializing EoS C library" << endl;
    initialize("../EoS_BQS_Derivatives/Coefficients_Parameters.dat");

    // load EoS tables, assess grid range
    std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    init_grid_ranges_only(quantityFile, derivFile);
  }
  else  // if thermo not from static C library, read in table from file
  {
    // initialize things needed to store eos table from file
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
    if ( use_nonconformal_extension )
    {
      tbqs_minima = {0.0,     -TBQS_INFINITY, -TBQS_INFINITY, -TBQS_INFINITY};
      tbqs_maxima = {TBQS_INFINITY, TBQS_INFINITY,  TBQS_INFINITY,  TBQS_INFINITY};
    }
    else
    {
      tbqs_minima = equation_of_state_table.get_grid_minima();
      tbqs_maxima = equation_of_state_table.get_grid_maxima();
    }

    // needed when using non-conformal extension to know actual table limits
    tbqs_minima_no_ext = equation_of_state_table.get_grid_minima();
    tbqs_maxima_no_ext = equation_of_state_table.get_grid_maxima();
  }

    // sets EoS type
    name = "table";

}





void EoS_table::init_grid_ranges_only(string quantityFile, string derivFile)
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

  // initialize grid ranges here
  if ( use_nonconformal_extension )
  {
    tbqs_minima = {0.0,          -TBQS_INFINITY, -TBQS_INFINITY, -TBQS_INFINITY};
    tbqs_maxima = {TBQS_INFINITY, TBQS_INFINITY,  TBQS_INFINITY,  TBQS_INFINITY};
  }
  else
  {
    tbqs_minima = {minT, minMuB, minMuQ, minMuS};
    tbqs_maxima = {maxT, maxMuB, maxMuQ, maxMuS};
  }

  // needed when using non-conformal extension to know actual table limits
  tbqs_minima_no_ext = {minT, minMuB, minMuQ, minMuS};
  tbqs_maxima_no_ext = {maxT, maxMuB, maxMuQ, maxMuS};

  dataFile.close();

	std::cout << "All initializations finished!" << std::endl;

  return;
}


////////////////////////////////////////////////////////////////////////////////
// point     = (T,muB,muQ,muS)
// densities = (e,rhoB,rhoS,rhoQ)
void EoS_table::get_eBSQ_densities_from_interpolator(
        double point[], double densities[] )  // point and densities both length = 4
{
    vector<double> results;
    equation_of_state_table.evaluate(
      vector<double>(point, point + 4), results,
      vector<string>({ "e","B","S","Q" }) );
    std::copy(results.begin(), results.end(), densities);
}

////////////////////////////////////////////////////////////////////////////////
// point     = (T,muB,muQ,muS)
// densities = (s,rhoB,rhoS,rhoQ)
void EoS_table::get_sBSQ_densities_from_interpolator(
        double point[], double densities[] )  // point and densities both length = 4
{
    vector<double> results;
    equation_of_state_table.evaluate(
      vector<double>(point, point + 4), results,
      vector<string>({ "s","B","S","Q" }) );
    std::copy(results.begin(), results.end(), densities);
}






////////////////////////////////////////////////////////////////////////////////
void EoS_table::get_eBSQ( double point_in[], double results[] )
                { get_eBSQ_safe( point_in, results ); }

void EoS_table::get_eBSQ_safe( const double point_in[], double results[] )
{
  double point_projected[4];
  for ( int i = 0; i < 4; i++ ) point_projected[i] = point_in[i];

//cout << __PRETTY_FUNCTION__ << ": " << point_projected[0] << "   " << point_projected[1] << "   "
//     << point_projected[2] << "   " << point_projected[3] << endl;

  //============================================================================
  // decide this w.r.t. the tbqs ranges sans extension
  // --> needed to decide whether to use extension
  bool point_not_in_range = !point_is_in_range_no_ext( point_projected[0], point_projected[1],
                                                       point_projected[2], point_projected[3] );

  //============================================================================
  if ( use_nonconformal_extension and point_not_in_range )
  {
    /// NOTE: point gets reset!
    // project back toward origin until intersecting grid boundary
    eos_extension::project_to_boundary(
        point_projected, tbqs_minima_no_ext.data(), tbqs_maxima_no_ext.data() );
  }

  //============================================================================
  // evaluate the relevant grid point
  if (use_static_C_library)
    STANDARD_get_eBSQ_densities(point_projected, results);
  else  // using table itself
    get_eBSQ_densities_from_interpolator(point_projected, results);

  //============================================================================
  // project back to original point using non-conformal extension
  if ( use_nonconformal_extension and point_not_in_range )
    eos_extension::get_nonconformal_extension( point_in, point_projected, results, 1 );
}






////////////////////////////////////////////////////////////////////////////////
void EoS_table::get_sBSQ( double point_in[], double results[] )
                { get_sBSQ_safe( point_in, results ); }

void EoS_table::get_sBSQ_safe( const double point_in[], double results[] )
{
  double point_projected[4];
  for ( int i = 0; i < 4; i++ ) point_projected[i] = point_in[i];

//cout << __PRETTY_FUNCTION__ << ": " << point_projected[0] << "   " << point_projected[1] << "   "
//     << point_projected[2] << "   " << point_projected[3] << endl;

  //============================================================================
  // decide this w.r.t. the tbqs ranges sans extension
  // --> needed to decide whether to use extension
  bool point_not_in_range = !point_is_in_range_no_ext( point_projected[0], point_projected[1],
                                                       point_projected[2], point_projected[3] );

  //============================================================================
  if ( use_nonconformal_extension and point_not_in_range )
  {
    /// NOTE: point gets reset!
    // project back toward origin until intersecting grid boundary
    eos_extension::project_to_boundary(
        point_projected, tbqs_minima_no_ext.data(), tbqs_maxima_no_ext.data() );
  }

  //============================================================================
  // evaluate the relevant grid point
  if (use_static_C_library)
    STANDARD_get_sBSQ_densities(point_projected, results);
  else  // using table itself
    get_sBSQ_densities_from_interpolator(point_projected, results);

  //============================================================================
  // project back to original point using non-conformal extension
  if ( use_nonconformal_extension and point_not_in_range )
    eos_extension::get_nonconformal_extension( point_in, point_projected, results, 2 );
}





////////////////////////////////////////////////////////////////////////////////
void EoS_table::get_full_thermo( double point_in[], double results[] )
                { get_full_thermo_safe( point_in, results ); }

void EoS_table::get_full_thermo_safe( const double point_in[], double results[] )
{
  double point_projected[4];
  for ( int i = 0; i < 4; i++ ) point_projected[i] = point_in[i];

//cout << __PRETTY_FUNCTION__ << ": " << point_projected[0] << "   " << point_projected[1] << "   "
//     << point_projected[2] << "   " << point_projected[3] << endl;

  //============================================================================
  // decide this w.r.t. the tbqs ranges sans extension
  // --> needed to decide whether to use extension
  bool point_not_in_range = !point_is_in_range_no_ext( point_projected[0], point_projected[1],
                                                       point_projected[2], point_projected[3] );

  //============================================================================
  if ( use_nonconformal_extension and point_not_in_range )
  {
    /// NOTE: point gets reset!
    // project back toward origin until intersecting grid boundary
    eos_extension::project_to_boundary(
        point_projected, tbqs_minima_no_ext.data(), tbqs_maxima_no_ext.data() );

const double hc = constants::hbarc_MeVfm;

cout << "Original point: " << point_in[0]*hc << "   " << point_in[1]*hc << "   "
      << point_in[2]*hc << "   " << point_in[3]*hc << endl;
cout << "Projected point: " << point_projected[0]*hc << "   " << point_projected[1]*hc << "   "
      << point_projected[2]*hc << "   " << point_projected[3]*hc << endl;
  }

  //============================================================================
  // evaluate the relevant grid point
  if (use_static_C_library)
    STANDARD_get_full_thermo( point_projected, results );
  else  // using table itself
  {
    // copy C arrays to C++ vectors
    vector<double> v_point(point_projected, point_projected+4);
    vector<double> v_results(results, results+17);

    // evaluate EoS interpolator at current location (S and Q NOT SWAPPED)
    equation_of_state_table.evaluate( v_point, v_results ); 

    if ( v_results.size() != 17 )
    {
      cerr << "PROBLEM" << endl;
      exit(1);
    }

    // copy C++ vector of results back to C array
    std::copy(v_results.begin(), v_results.end(), results);
  }


  //============================================================================
  // project back to original point using non-conformal extension
  if ( use_nonconformal_extension and point_not_in_range )
    eos_extension::get_nonconformal_extension( point_in, point_projected, results, 0 );

//cout << "Thermo:" << endl;
//for (int i = 0; i < 17; i++)
//  cout << results[i] << endl;

}  


////////////////////////////////////////////////////////////////////////////////
bool EoS_table::point_is_in_range_no_ext(
                double setT, double setmuB, double setmuQ, double setmuS )
{
//cout << __PRETTY_FUNCTION__ << ": " << setT << "   " << setmuB << "   "
//     << setmuQ << "   " << setmuS << endl;
  if(setT < tbqs_minima_no_ext[0] || setT > tbqs_maxima_no_ext[0])
  {
    if ( VERBOSE > 3 )
    std::cout << "T = " << setT
      << " is out of table range (ignoring extension)."
         " Valid values are between ["
      << tbqs_minima_no_ext[0] << "," << tbqs_maxima_no_ext[0] << "]" << std::endl;
    return false;
  }
  if(setmuB < tbqs_minima_no_ext[1] || setmuB > tbqs_maxima_no_ext[1])
  {
    if ( VERBOSE > 3 )
    std::cout << "muB = " << setmuB
      << " is out of table range (ignoring extension)."
         " Valid values are between ["
      << tbqs_minima_no_ext[1] << "," << tbqs_maxima_no_ext[1] << "]" << std::endl;
    return false;
  }
  if(setmuQ < tbqs_minima_no_ext[2] || setmuQ > tbqs_maxima_no_ext[2])
  {
    if ( VERBOSE > 3 )
    std::cout << "muQ = " << setmuQ
      << " is out of table range (ignoring extension)."
         " Valid values are between ["
      << tbqs_minima_no_ext[2] << "," << tbqs_maxima_no_ext[2] << "]" << std::endl;
    return false;
  }
  if(setmuS < tbqs_minima_no_ext[3] || setmuS > tbqs_maxima_no_ext[3])
  {
    if ( VERBOSE > 3 )
    std::cout << "muS = " << setmuS
      << " is out of table range (ignoring extension)."
         " Valid values are between ["
      << tbqs_minima_no_ext[3] << "," << tbqs_maxima_no_ext[3] << "]" << std::endl;
    return false;
  }
  return true;
}