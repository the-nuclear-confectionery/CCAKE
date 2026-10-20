#include "../include/eos_table.h"

#include <algorithm>
#include <functional>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include "../include/eos_extension.h"
#include "../include/eos_header.h"
#include "../include/eos_table.h"

EoS_table::EoS_table( string eos_path, bool normalize_by_T )
{
  equation_of_state_table = std::make_shared<InterpolatorND<4>>();
  //////////////////////////////////////////////////////////////////////////////
  // read in table from file
  // initialize things needed to store eos table from file
  if ( using_HDF )
    equation_of_state_table_filename = eos_path + "/thermo.h5";
  else
    equation_of_state_table_filename = eos_path + "/thermo.dat";

  cout << "Initializing EoS from input file(s): "
       << equation_of_state_table_filename << endl;
  equation_of_state_table->initialize( equation_of_state_table_filename );

  // set names of EoS quantities to interpolate, in order
  equation_of_state_table->set_grid_names(
    vector<string>{ "T","muB","muQ","muS" } );
  // Auto-detect the optional 18th field chiBB_nc: the non-critical chi_BB baseline, a smooth
  // fit to chi_BB of THIS table with the critical ridge masked out.  When present it is used as
  // xi_0 in critical_scaling_factor_bulk instead of the global constant, so the enhancement is
  // measured against the local non-critical value.  Tables written before this (21 columns) are
  // read exactly as before.
  {
    // The file opens with a 4-column grid header, so probing only the first non-comment line
    // reports 4 and silently selects the 21-column layout.  Take the widest of the first few
    // lines instead, which is a data row in every table we ship.
    std::ifstream probe( equation_of_state_table_filename );
    std::string line;
    int ncol = 0, examined = 0;
    while ( examined < 6 && std::getline( probe, line ) )
    {
      if ( line.empty() || line[0] == '#' ) continue;
      ++examined;
      std::istringstream iss( line );
      int n = 0; double dummy;
      while ( iss >> dummy ) ++n;
      if ( n > ncol ) ncol = n;
    }
    has_chiBB_nc = ( ncol >= 22 );
    cout << "EoS table has " << ncol << " columns; chiBB_nc "
         << ( has_chiBB_nc ? "PRESENT (local xi_0)" : "absent (global xi_0)" ) << endl;
  }

  vector<string> field_names{ "p","s","B","S","Q","e","cs2",
                              "chiBB","chiQQ","chiSS","chiBQ","chiBS",
                              "chiQS","chiTB","chiTQ","chiTS","chiTT" };
  if ( has_chiBB_nc ) field_names.push_back( "chiBB_nc" );
  equation_of_state_table->set_field_names( field_names );

  // finally, all dimensionful quantities should be convert to fm and
  // all dimensionless quantities need to be rescaled appropriately
  /*equation_of_state_table.rescale_axis( "T",   1.0/hbarc_MeVfm );
  equation_of_state_table.rescale_axis( "muB", 1.0/hbarc_MeVfm );
  equation_of_state_table.rescale_axis( "muQ", 1.0/hbarc_MeVfm );
  equation_of_state_table.rescale_axis( "muS", 1.0/hbarc_MeVfm );*/
  equation_of_state_table->rescale_axes( 1.0/hbarc_MeVfm );  // do all axes at once

  if ( normalize_by_T )
  {
    // Table stores dimensionless ratios (p/T^4, s/T^3, etc.);
    // multiply by the appropriate power of T (already in fm^-1) to get physical units.
    equation_of_state_table->rescale( "p",     "T", 4 );
    equation_of_state_table->rescale( "e",     "T", 4 );
    equation_of_state_table->rescale( "s",     "T", 3 );
    equation_of_state_table->rescale( "B",     "T", 3 );
    equation_of_state_table->rescale( "S",     "T", 3 );
    equation_of_state_table->rescale( "Q",     "T", 3 );
    equation_of_state_table->rescale( "chiBB", "T", 2 );
    if ( has_chiBB_nc ) equation_of_state_table->rescale( "chiBB_nc", "T", 2 );
    equation_of_state_table->rescale( "chiQQ", "T", 2 );
    equation_of_state_table->rescale( "chiSS", "T", 2 );
    equation_of_state_table->rescale( "chiBQ", "T", 2 );
    equation_of_state_table->rescale( "chiBS", "T", 2 );
    equation_of_state_table->rescale( "chiQS", "T", 2 );
    equation_of_state_table->rescale( "chiTB", "T", 2 );
    equation_of_state_table->rescale( "chiTQ", "T", 2 );
    equation_of_state_table->rescale( "chiTS", "T", 2 );
    equation_of_state_table->rescale( "chiTT", "T", 2 );
  }
  else
  {
    // Table stores fields in MeV^n; multiply by (1/hbarc_MeVfm)^n to convert to fm^-n.
    const double hc1 = 1.0/hbarc_MeVfm;
    const double hc2 = hc1*hc1;
    const double hc3 = hc2*hc1;
    const double hc4 = hc3*hc1;
    equation_of_state_table->rescale_field( "p",     hc4 );
    equation_of_state_table->rescale_field( "e",     hc4 );
    equation_of_state_table->rescale_field( "s",     hc3 );
    equation_of_state_table->rescale_field( "B",     hc3 );
    equation_of_state_table->rescale_field( "S",     hc3 );
    equation_of_state_table->rescale_field( "Q",     hc3 );
    equation_of_state_table->rescale_field( "chiBB", hc2 );
    if ( has_chiBB_nc ) equation_of_state_table->rescale_field( "chiBB_nc", hc2 );
    equation_of_state_table->rescale_field( "chiQQ", hc2 );
    equation_of_state_table->rescale_field( "chiSS", hc2 );
    equation_of_state_table->rescale_field( "chiBQ", hc2 );
    equation_of_state_table->rescale_field( "chiBS", hc2 );
    equation_of_state_table->rescale_field( "chiQS", hc2 );
    equation_of_state_table->rescale_field( "chiTB", hc2 );
    equation_of_state_table->rescale_field( "chiTQ", hc2 );
    equation_of_state_table->rescale_field( "chiTS", hc2 );
    equation_of_state_table->rescale_field( "chiTT", hc2 );
  }



  // set grid ranges
  if ( use_nonconformal_extension )
  {
    tbqs_minima = {0.0,     -TBQS_INFINITY, -TBQS_INFINITY, -TBQS_INFINITY};
    tbqs_maxima = {TBQS_INFINITY, TBQS_INFINITY,  TBQS_INFINITY,  TBQS_INFINITY};
  }
  else
  {
    tbqs_minima = equation_of_state_table->get_grid_minima();
    tbqs_maxima = equation_of_state_table->get_grid_maxima();
  }

  // needed when using non-conformal extension to know actual table limits
  tbqs_minima_no_ext = equation_of_state_table->get_grid_minima();
  tbqs_maxima_no_ext = equation_of_state_table->get_grid_maxima();

  //print the tbqs ranges
  cout << "T range: " << tbqs_minima[0] << "   " << tbqs_maxima[0] << endl;
  cout << "muB range: " << tbqs_minima[1] << "   " << tbqs_maxima[1] << endl;
  cout << "muS range: " << tbqs_minima[2] << "   " << tbqs_maxima[2] << endl;
  cout << "muQ range: " << tbqs_minima[3] << "   " << tbqs_maxima[3] << endl;


    // sets EoS type
    name = "table";

}




////////////////////////////////////////////////////////////////////////////////
// point     = (T,muB,muQ,muS)
// densities = (e,rhoB,rhoS,rhoQ)
void EoS_table::get_eBSQ_densities_from_interpolator(
        double point[], double densities[] )  // point and densities both length = 4
{
    static thread_local vector<double> tl_point(4), tl_results;
    static const vector<string> fields({ "e","B","S","Q" });
    std::copy(point, point + 4, tl_point.begin());
    equation_of_state_table->evaluate( tl_point, tl_results, fields );
    std::copy(tl_results.begin(), tl_results.end(), densities);
}

////////////////////////////////////////////////////////////////////////////////
// point     = (T,muB,muQ,muS)
// densities = (s,rhoB,rhoS,rhoQ)
void EoS_table::get_sBSQ_densities_from_interpolator(
        double point[], double densities[] )  // point and densities both length = 4
{
    static thread_local vector<double> tl_point(4), tl_results;
    static const vector<string> fields({ "s","B","S","Q" });
    std::copy(point, point + 4, tl_point.begin());
    equation_of_state_table->evaluate( tl_point, tl_results, fields );
    std::copy(tl_results.begin(), tl_results.end(), densities);
}






////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void EoS_table::get_eBSQ( double point_in[], double results[] )
                { get_eBSQ_safe( point_in, results ); }

void EoS_table::get_eBSQ_safe( const double point_in[], double results[] )
{
  //============================================================================
  // save input point in case we need to project to boundary
  // full results needed if projecting with extension
  double point_projected[4], results_full[18];
  for ( int i = 0; i < 4; i++ ) point_projected[i] = point_in[i];


  //============================================================================
  // decide this w.r.t. the tbqs ranges sans extension
  // --> needed to decide whether to use extension
  bool point_not_in_range
        = !point_is_in_range_no_ext( point_projected[0], point_projected[1],
                                     point_projected[2], point_projected[3] );


  //============================================================================
  if ( use_nonconformal_extension and point_not_in_range )
  {
    /// NOTE: point gets reset!
    // project back toward origin until intersecting grid boundary
    eos_extension::project_to_boundary(
        point_projected, tbqs_minima_no_ext.data(), tbqs_maxima_no_ext.data() );

    if ( VERBOSE > 8 )
    {
      const double hc = constants::hbarc_MeVfm;
      std::cout << "Original point: "
                << point_in[0]*hc << "   "
                << point_in[1]*hc << "   "
                << point_in[2]*hc << "   "
                << point_in[3]*hc << std::endl;
      std::cout << "Projected point: "
                << point_projected[0]*hc << "   "
                << point_projected[1]*hc << "   "
                << point_projected[2]*hc << "   "
                << point_projected[3]*hc << std::endl;
    }

    //============================================================================
    // MUST USE FULL THERMO TO SET NON-CONFORMAL EXTENSION
    // evaluate the relevant grid point
    // using table itself
    {
      static thread_local vector<double> v_point(4), v_results;
      std::copy(point_projected, point_projected+4, v_point.begin());
      equation_of_state_table->evaluate( v_point, v_results );

      // 18 when the table carries the optional chiBB_nc baseline, 17 otherwise
      const size_t n_expected = has_chiBB_nc ? 18 : 17;
      if ( v_results.size() != n_expected )
      {
        cerr << "PROBLEM: EoS table returned " << v_results.size()
             << " fields, expected " << n_expected << endl;
        exit(1);
      }
      std::copy(v_results.begin(), v_results.end(), results_full);
    }

    // project back to original point using non-conformal extension
    eos_extension::get_nonconformal_extension( point_in, point_projected,
                                               results_full );

    // option to print results
    if ( VERBOSE > 8 )
    {
      cout << "Thermo:" << endl;
      for (int i = 0; i < 17; i++)
        cout << results_full[i] << endl;
    }

    // set relevant densities from full thermo and return
    results[0] = results_full[5];
    results[1] = results_full[2];
    results[2] = results_full[3];
    results[3] = results_full[4];

  }
  else
  {
    //============================================================================
    // evaluate the relevant grid point
    // using table itself
    get_eBSQ_densities_from_interpolator(point_projected, results);
  }

}






////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void EoS_table::get_sBSQ( double point_in[], double results[] )
                { get_sBSQ_safe( point_in, results ); }

void EoS_table::get_sBSQ_safe( const double point_in[], double results[] )
{
  //============================================================================
  // save input point in case we need to project to boundary
  // full results needed if projecting with extension
  double point_projected[4], results_full[18];
  for ( int i = 0; i < 4; i++ ) point_projected[i] = point_in[i];


  //============================================================================
  // decide this w.r.t. the tbqs ranges sans extension
  // --> needed to decide whether to use extension
  bool point_not_in_range
        = !point_is_in_range_no_ext( point_projected[0], point_projected[1],
                                     point_projected[2], point_projected[3] );


  //============================================================================
  if ( use_nonconformal_extension and point_not_in_range )
  {
    /// NOTE: point gets reset!
    // project back toward origin until intersecting grid boundary
    eos_extension::project_to_boundary(
        point_projected, tbqs_minima_no_ext.data(), tbqs_maxima_no_ext.data() );

    if ( VERBOSE > 8 )
    {
      const double hc = constants::hbarc_MeVfm;
      std::cout << "Original point: "
                << point_in[0]*hc << "   "
                << point_in[1]*hc << "   "
                << point_in[2]*hc << "   "
                << point_in[3]*hc << std::endl;
      std::cout << "Projected point: "
                << point_projected[0]*hc << "   "
                << point_projected[1]*hc << "   "
                << point_projected[2]*hc << "   "
                << point_projected[3]*hc << std::endl;
    }

    //============================================================================
    // MUST USE FULL THERMO TO SET NON-CONFORMAL EXTENSION
    // evaluate the relevant grid point
    // using table itself
    {
      static thread_local vector<double> v_point(4), v_results;
      std::copy(point_projected, point_projected+4, v_point.begin());
      equation_of_state_table->evaluate( v_point, v_results );

      // 18 when the table carries the optional chiBB_nc baseline, 17 otherwise
      const size_t n_expected = has_chiBB_nc ? 18 : 17;
      if ( v_results.size() != n_expected )
      {
        cerr << "PROBLEM: EoS table returned " << v_results.size()
             << " fields, expected " << n_expected << endl;
        exit(1);
      }
      std::copy(v_results.begin(), v_results.end(), results_full);
    }

    // project back to original point using non-conformal extension
    eos_extension::get_nonconformal_extension( point_in, point_projected,
                                               results_full );

    // option to print results
    if ( VERBOSE > 8 )
    {
      cout << "Thermo:" << endl;
      for (int i = 0; i < 17; i++)
        cout << results_full[i] << endl;
    }

    // set relevant densities from full thermo and return
    results[0] = results_full[1];
    results[1] = results_full[2];
    results[2] = results_full[3];
    results[3] = results_full[4];

  }
  else
  {
    //============================================================================
    // evaluate the relevant grid point
    // using table itself
    get_sBSQ_densities_from_interpolator(point_projected, results);
  }

}





////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void EoS_table::get_full_thermo( double point_in[], double results[] )
                { get_full_thermo_safe( point_in, results ); }

void EoS_table::get_full_thermo_safe( const double point_in[], double results[] )
{
  //============================================================================
  // save input point in case we need to project to boundary
  double point_projected[4];
  for ( int i = 0; i < 4; i++ ) point_projected[i] = point_in[i];

  //============================================================================
  // decide this w.r.t. the tbqs ranges sans extension
  // --> needed to decide whether to use extension
  bool point_not_in_range
        = !point_is_in_range_no_ext( point_projected[0], point_projected[1],
                                     point_projected[2], point_projected[3] );

  //============================================================================
  if ( use_nonconformal_extension and point_not_in_range )
  {
    /// NOTE: point gets reset!
    // project back toward origin until intersecting grid boundary
    eos_extension::project_to_boundary(
        point_projected, tbqs_minima_no_ext.data(), tbqs_maxima_no_ext.data() );

    // print if needed
    if ( VERBOSE > 8 )
    {
      const double hc = constants::hbarc_MeVfm;
      std::cout << "Original point: "
                << point_in[0]*hc << "   "
                << point_in[1]*hc << "   "
                << point_in[2]*hc << "   "
                << point_in[3]*hc << std::endl;
      std::cout << "Projected point: "
                << point_projected[0]*hc << "   "
                << point_projected[1]*hc << "   "
                << point_projected[2]*hc << "   "
                << point_projected[3]*hc << std::endl;
    }
  }

  //============================================================================
  // evaluate the relevant grid point
  // using table itself
  {
    static thread_local vector<double> v_point(4), v_results;
    std::copy(point_projected, point_projected+4, v_point.begin());
    equation_of_state_table->evaluate( v_point, v_results );

    // 18 when the table carries the optional chiBB_nc baseline, 17 otherwise
    const size_t n_expected = has_chiBB_nc ? 18 : 17;
    if ( v_results.size() != n_expected )
    {
      cerr << "PROBLEM: EoS table returned " << v_results.size()
           << " fields, expected " << n_expected << endl;
      exit(1);
    }
    std::copy(v_results.begin(), v_results.end(), results);
  }


  //============================================================================
  // project back to original point using non-conformal extension
  if ( use_nonconformal_extension and point_not_in_range )
    eos_extension::get_nonconformal_extension( point_in, point_projected,
                                               results );
  // get_nonconformal_extension rewrites fields 0..16 but leaves field 17 holding chiBB_nc at
  // the BOUNDARY point, which no longer pairs with the extension's own chiBB.  Set the local
  // baseline equal to chiBB there so the ratio is exactly 1: no critical enhancement outside
  // the table, which is where the critical point is not.
  if ( has_chiBB_nc and use_nonconformal_extension and point_not_in_range )
    results[17] = results[7];

  // option to print results
  if ( VERBOSE > 8 )
  {
    cout << "Thermo:" << " ";
    for (int i = 0; i < 17; i++)
      cout << results[i] << " ";
    cout << endl;
  }

}  


////////////////////////////////////////////////////////////////////////////////
bool EoS_table::point_is_in_range_no_ext(
                double setT, double setmuB, double setmuQ, double setmuS )
{
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