#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
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
#include "../src/eos_derivatives.cpp"
#include "../src/eos_initialization.cpp"

using namespace constants;
using namespace ccake;
using std::vector;
using std::string;

////////////////////////////////////////////////////////////////////////////////
// Pre-allocate GSL workspace for derivative computations
EquationOfState::EquationOfState() {
  _dv    = gsl_vector_alloc(3);
  _da    = gsl_vector_alloc(3);
  _db    = gsl_vector_alloc(3);
  _dm    = gsl_matrix_alloc(3, 3);
  _dp    = gsl_permutation_alloc(3);
  _dminv = gsl_matrix_alloc(3, 3);
  _dy    = gsl_vector_alloc(3);
}
EquationOfState::~EquationOfState() {
  if (_dv)    gsl_vector_free(_dv);
  if (_da)    gsl_vector_free(_da);
  if (_db)    gsl_vector_free(_db);
  if (_dm)    gsl_matrix_free(_dm);
  if (_dp)    gsl_permutation_free(_dp);
  if (_dminv) gsl_matrix_free(_dminv);
  if (_dy)    gsl_vector_free(_dy);
}

////////////////////////////////////////////////////////////////////////////////

vector<double> EquationOfState::get_thermodynamics( vector<double> & tbqsIn,
                                                      const string & eos_name )
{
  if ( eos_name == "default" )
    tbqs( tbqsIn, chosen_EOS_map[default_eos_name] );
  else
    tbqs( tbqsIn, chosen_EOS_map[eos_name] );

  return vector<double>({ pVal,entrVal,BVal,SVal,QVal,eVal,cs2Val,
                          db2,dq2,ds2,dbdq,dbds,dsdq,dtdb,dtdq,dtds,dt2 });
}



////////////////////////////////////////////////////////////////////////////////
bool EquationOfState::point_not_in_range(
                        double setT, double setmuB, double setmuQ,
                        double setmuS, pEoS_base peos )
{
  const auto & tbqs_mins = peos->tbqs_minima;
  const auto & tbqs_maxs = peos->tbqs_maxima;

  if(setT < tbqs_mins[0] || setT > tbqs_maxs[0])
  { 
    std::cout << "T = " << setT
      << " is out of (" << peos->name << ") range."
         " Valid values are between ["
      << tbqs_mins[0] << "," << tbqs_maxs[0] << "]" << std::endl;
    return true;
  }
  if(setmuB < tbqs_mins[1] || setmuB > tbqs_maxs[1])
  {
    std::cout << "muB = " << setmuB
      << " is out of (" << peos->name << ") range."
         " Valid values are between ["
      << tbqs_mins[1] << "," << tbqs_maxs[1] << "]" << std::endl;
    return true;
  }
  if(setmuQ < tbqs_mins[2] || setmuQ > tbqs_maxs[2])
  {
    std::cout << "muQ = " << setmuQ
      << " is out of (" << peos->name << ") range."
         " Valid values are between ["
      << tbqs_mins[2] << "," << tbqs_maxs[2] << "]" << std::endl;
    return true;
  }
  if(setmuS < tbqs_mins[3] || setmuS > tbqs_maxs[3])
  {
    std::cout << "muS = " << setmuS
      << " is out of (" << peos->name << ") range."
         " Valid values are between ["
      << tbqs_mins[3] << "," << tbqs_maxs[3] << "]" << std::endl;
    return true;
  }
  return false;
}




void EquationOfState::tbqs( double setT, double setmuB, double setmuQ, double setmuS,
                            const string & eos_name )
{
  if ( eos_name == "default" )
    tbqs( setT, setmuB, setmuQ, setmuS, chosen_EOS_map[default_eos_name] );
  else
    tbqs( setT, setmuB, setmuQ, setmuS, chosen_EOS_map[eos_name] );
}

void EquationOfState::tbqs( vector<double> & tbqsIn, const string & eos_name )
{
  if ( eos_name == "default" )
    tbqs( tbqsIn, chosen_EOS_map[default_eos_name] );
  else
    tbqs( tbqsIn, chosen_EOS_map[eos_name] );
}



////////////////////////////////////////////////////////////////////////////////
void EquationOfState::tbqs( double setT, double setmuB, double setmuQ,
                            double setmuS, pEoS_base peos )
{
  // set name of current EoS
  current_eos_name = peos->name;
  // Detect degenerate axes from the EoS coverage (min==max means "frozen")
  auto frozen = [](double a, double b) { return std::abs(a - b) < 1e-30; };

  const bool muB_frozen = frozen(peos->tbqs_minima[1], peos->tbqs_maxima[1]);
  const bool muQ_frozen = frozen(peos->tbqs_minima[2], peos->tbqs_maxima[2]);
  const bool muS_frozen = frozen(peos->tbqs_minima[3], peos->tbqs_maxima[3]);

  T_only_mode      = (muB_frozen && muQ_frozen && muS_frozen);
  baryon_only_mode = (!T_only_mode && muQ_frozen && muS_frozen);
  // check if proposed point is in range
  bool point_is_in_range = !point_not_in_range( setT, setmuB, setmuQ, setmuS, peos );
  if ( point_is_in_range )
  {
    tbqsPosition[0] = setT;
    tbqsPosition[1] = setmuB;
    tbqsPosition[2] = setmuQ;
    tbqsPosition[3] = setmuS;

    // if we are in range, compute all thermodynamic quantities at this point
    evaluate_thermodynamics( peos );
  }
  else
  {
    std::cout << "EoS: " << peos->name << " does not cover this point. You need to tell me what to do!" << std::endl;
    std::cerr << "Point was out of range! You need to tell me what to do!" << std::endl;
    exit(8);
  }

  return;
}


void EquationOfState::evaluate_thermodynamics( pEoS_base peos )
{
  //============================================================================
  // this function now works the same for all EoSs by construction
  vector<double> thermodynamics;
  double phase_diagram_point[4] = { tbqsPosition[0], tbqsPosition[1],
                                    tbqsPosition[2], tbqsPosition[3] };

  double thermo_array[17];
  peos->get_full_thermo( phase_diagram_point, thermo_array );

  thermodynamics.assign(thermo_array, thermo_array + 17);

  //============================================================================
  // set final thermodynamic results
  pVal    = thermodynamics[0];
  entrVal = thermodynamics[1];
  BVal    = thermodynamics[2];
  SVal    = thermodynamics[3];
  QVal    = thermodynamics[4];
  eVal    = thermodynamics[5];
  cs2Val  = thermodynamics[6];
  db2     = thermodynamics[7];
  dq2     = thermodynamics[8];
  ds2     = thermodynamics[9];
  dbdq    = thermodynamics[10];
  dbds    = thermodynamics[11];
  dsdq    = thermodynamics[12];
  dtdb    = thermodynamics[13];
  dtdq    = thermodynamics[14];
  dtds    = thermodynamics[15];
  dt2     = thermodynamics[16];


  //recompute the EoS but assuming muQ = muB = muS = 0 to get chiBB0
  if ( muB() != 0.0 || muQ() != 0.0 || muS() != 0.0 )
  {    double zero_mu_phase_diagram_point[4] = { tbqsPosition[0], 0.0, 0.0, 0.0 };
       double zero_mu_thermo_array[17];
       peos->get_full_thermo( zero_mu_phase_diagram_point, zero_mu_thermo_array );
       db20 = zero_mu_thermo_array[7];
  }
  else  {
    db20 = db2;
  }

  //if EoS is not table, set chiBB0 = chiBB
  if ( peos->name != "table" ) db20 = db2;
//cout << "before" << endl;
//cout << "THERMO DUMP: " << pVal << "   " << entrVal << "   " << BVal << "   "
//      << SVal << "   " << QVal << "   " << eVal << "   " << cs2Val << "   "
//      << db2 << "   " << dq2 << "   " << ds2 << "   " << dbdq << "   "
//      << dbds << "   " << dsdq << "   " << dtdb << "   " << dtdq << "   "
//      << dtds << "   " << dt2 << "   " << peos->name << endl;
//cout << "after" << endl;

/*if (true) exit(1);*/
}



////////////////////////////////////////////////////////////////////////////////
// Set cached thermodynamics directly from a gap-file row.
// current_eos_name is set to "table" so the skip-list in the next timestep's
// rootfinder loop still starts at the table EoS rather than conformal fallbacks.
void EquationOfState::tbqs_from_gap(double T_g, double muB_g, double muQ_g, double muS_g,
                                     const double thermo_g[N_GAP_THERMO])
{
    tbqsPosition[0] = T_g;
    tbqsPosition[1] = muB_g;
    tbqsPosition[2] = muQ_g;
    tbqsPosition[3] = muS_g;

    pVal    = thermo_g[0];
    entrVal = thermo_g[1];
    BVal    = thermo_g[2];
    SVal    = thermo_g[3];
    QVal    = thermo_g[4];
    eVal    = thermo_g[5];
    cs2Val  = thermo_g[6];
    db2     = thermo_g[7];
    dq2     = thermo_g[8];
    ds2     = thermo_g[9];
    dbdq    = thermo_g[10];
    dbds    = thermo_g[11];
    dsdq    = thermo_g[12];
    dtdb    = thermo_g[13];
    dtdq    = thermo_g[14];
    dtds    = thermo_g[15];
    dt2     = thermo_g[16];
    db20    = db2;  // gap table doesn't carry the mu=0 susceptibility separately

    current_eos_name = "table";
}

double EquationOfState::T()   const { return tbqsPosition[0]; }
double EquationOfState::muB() const { return tbqsPosition[1]; }
double EquationOfState::muQ() const { return tbqsPosition[2]; }
double EquationOfState::muS() const { return tbqsPosition[3]; }

double EquationOfState::p()   const { return pVal; }
double EquationOfState::s()   const { return entrVal; }
double EquationOfState::B()   const { return BVal; }
double EquationOfState::S()   const { return SVal; }
double EquationOfState::Q()   const { return QVal; }
double EquationOfState::e()   const { return eVal; }
double EquationOfState::cs2() const { return cs2Val; }
double EquationOfState::w()   const { return eVal + pVal; }

double EquationOfState::chiBB() const { return db2; }
double EquationOfState::chiBB0() const { return db20; }

double EquationOfState::dwds()
{
  // include non-zero charge terms
	double charge_terms	= 0.0;
  if ( abs(BVal)>TINY ) charge_terms += BVal/dentr_dmub();
  if ( abs(QVal)>TINY ) charge_terms += QVal/dentr_dmuq();
  if ( abs(SVal)>TINY ) charge_terms += SVal/dentr_dmus();

	/*if ( check_derivatives )
	cout << endl << endl << "inside dwds(): "
		<< T() << "   " << entrVal << "   " << dentr_dt() << "   "
		<< BVal << "   " << dentr_dmub() << "   "
		<< SVal << "   " << dentr_dmus() << "   "
		<< QVal << "   " << dentr_dmuq() << endl << endl;*/

  return T() + entrVal/dentr_dt() + charge_terms;
}

double EquationOfState::dwdB()
{
  // include non-zero charge terms
	double charge_terms	= 0.0;
  if ( abs(BVal)>TINY )
  {
    charge_terms += entrVal/db_dt() + BVal/db_dmub();
    if ( abs(QVal)>TINY ) charge_terms += QVal/db_dmuq();
    if ( abs(SVal)>TINY ) charge_terms += SVal/db_dmus();
  }

  return muB() + charge_terms;
}

double EquationOfState::dwdS()
{
  // include non-zero charge terms
	double charge_terms	= 0.0;
  if ( abs(SVal)>TINY )
  {
    charge_terms += entrVal/ds_dt() + SVal/ds_dmus();
    if ( abs(BVal)>TINY ) charge_terms += BVal/ds_dmub();
    if ( abs(QVal)>TINY ) charge_terms += QVal/ds_dmuq();
  }

  return muS() + charge_terms;
}

double EquationOfState::dwdQ()
{
  // include non-zero charge terms
	double charge_terms	= 0.0;
  if ( abs(QVal)>TINY )
  {
    charge_terms += entrVal/dq_dt() + QVal/dq_dmuq();
    if ( abs(BVal)>TINY ) charge_terms += BVal/dq_dmub();
    if ( abs(SVal)>TINY ) charge_terms += SVal/dq_dmus();
  }

  return muQ() + charge_terms;
}


double EquationOfState::dalpha_Bds()
{
  double charge_terms	= 0.0;
  if ( abs(BVal)>TINY )
  {
    charge_terms += 1/T()/dentr_dmub();
  }
  return  -muB()/T()/T()/dentr_dt() + charge_terms;
}

double EquationOfState::dalpha_Sds()
{
  double charge_terms	= 0.0;
  if ( abs(SVal)>TINY )
  {
    charge_terms += 1/T()/dentr_dmus();
  }

  return -muS()/T()/T()/dentr_dt() + charge_terms;
}

double EquationOfState::dalpha_Qds()
{
  double charge_terms	= 0.0;
  if ( abs(QVal)>TINY )
  {
    charge_terms += 1/T()/dentr_dmuq();
  }

  return -muQ()/T()/T()/dentr_dt() + charge_terms;
}


double EquationOfState::dalpha_BdB()
{
  double charge_terms	= 0.0;
  if ( abs(BVal)>TINY )
  {
    charge_terms += 1/T()/db_dmub() -muB()/T()/T()/db_dt();
  }
  return charge_terms;
}

double EquationOfState::dalpha_SdB()
{
  double charge_terms	= 0.0;
  if ( abs(BVal)>TINY )
  {
    charge_terms += -muS()/T()/T()/db_dt();
    if ( abs(SVal)>TINY )
    {
      charge_terms += 1/T()/db_dmus(); 
    }
  }


  return charge_terms;
}

double EquationOfState::dalpha_QdB()
{
  double charge_terms	= 0.0;
  if ( abs(BVal)>TINY )
  {
    charge_terms += -muQ()/T()/T()/db_dt();
    if ( abs(QVal)>TINY )
    {
      charge_terms +=  1/T()/db_dmuq();
    }
  }
  return charge_terms;
}

double EquationOfState::dalpha_BdS()
{
  double charge_terms	= 0.0;
  if ( abs(SVal)>TINY )
  {
    charge_terms += -muB()/T()/T()/ds_dt();
    if ( abs(BVal)>TINY )
    {
      charge_terms += 1/T()/ds_dmub();
    }
  }
  return charge_terms;
}

double EquationOfState::dalpha_SdS()
{
  double charge_terms	= 0.0;
  if ( abs(SVal)>TINY )
  {
    charge_terms += -muS()/T()/T()/ds_dt();
    if ( abs(SVal)>TINY )
    {
      charge_terms += 1/T()/ds_dmus();
    }
  }
  return charge_terms;
}

double EquationOfState::dalpha_QdS()
{
  double charge_terms	= 0.0;
  if ( abs(SVal)>TINY )
  {
    charge_terms += -muQ()/T()/T()/ds_dt();
    if ( abs(QVal)>TINY )
    {
      charge_terms += 1/T()/ds_dmuq();
    }
  }
  return charge_terms;
}

double EquationOfState::dalpha_BdQ()
{
  double charge_terms	= 0.0;
  if ( abs(QVal)>TINY )
  {
    charge_terms += -muB()/T()/T()/dq_dt();
    if ( abs(BVal)>TINY )
    {
      charge_terms += 1/T()/dq_dmub();
    }
  }
  return charge_terms;
}

double EquationOfState::dalpha_SdQ()
{
  double charge_terms	= 0.0;
  if ( abs(QVal)>TINY )
  {
    charge_terms += -muS()/T()/T()/dq_dt();
    if ( abs(SVal)>TINY )
    {
      charge_terms += 1/T()/dq_dmus();
    }
  }
  return charge_terms;
}

double EquationOfState::dalpha_QdQ()
{
  double charge_terms	= 0.0;
  if ( abs(QVal)>TINY )
  {
    charge_terms += -muQ()/T()/T()/dq_dt();
    if ( abs(QVal)>TINY )
    {
      charge_terms += 1/T()/dq_dmuq();
    }
  }
  return charge_terms;
}


double EquationOfState::cs2out(double Tt, const string & eos_name)
{  //return cs2 given t and mu's=0
//  std::cout << __PRETTY_FUNCTION__ << "::" << __LINE__ << ": "
//            << "using " << eos_name << "!" << std::endl;
  tbqs(Tt, 0.0, 0.0, 0.0, eos_name);
  return cs2Val;
}

double EquationOfState::cs2out(double Tt, double muBin, double muQin, double muSin, const string & eos_name)
{  //return cs2 given t and mu's
//  std::cout << __PRETTY_FUNCTION__ << "::" << __LINE__ << ": "
//            << "using " << eos_name << "!" << std::endl;
  tbqs(Tt, muBin, muQin, muSin, eos_name);
  return cs2Val;
}

double EquationOfState::wfz(double Tt, const string & eos_name)
{   // return e + p for tbqs
//  std::cout << __PRETTY_FUNCTION__ << "::" << __LINE__ << ": "
//            << "using " << eos_name << "!" << std::endl;
  tbqs(Tt, 0.0, 0.0, 0.0, eos_name);
  return eVal + pVal;
}

double EquationOfState::wfz(double Tt, double muBin, double muQin, double muSin, const string & eos_name)
{   // return e + p for tbqs
//  std::cout << __PRETTY_FUNCTION__ << "::" << __LINE__ << ": "
//            << "using " << eos_name << "!" << std::endl;
  tbqs(Tt, muBin, muQin, muSin, eos_name);
  return eVal + pVal;
}


double EquationOfState::s_terms_T(double Tt, const string & eos_name)
{
//  std::cout << __PRETTY_FUNCTION__ << "::" << __LINE__ << ": "
//            << "using " << eos_name << "!" << std::endl;
  tbqs(Tt, 0, 0, 0, eos_name);
  return entrVal;
}


double EquationOfState::s_terms_T( double Tt, double muBin,
                                   double muQin, double muSin,
                                   const string & eos_name )
{
  tbqs(Tt, muBin, muQin, muSin, eos_name);
  return entrVal;
}


// UNCOMMENTED BY C. PLUMBERG
void EquationOfState::eosin(std::string type){}

double EquationOfState::A() { return w()-s()*dwds(); }


// confirm with Jaki
double EquationOfState::efreeze(double T_freeze_out_at_mu_eq_0, const string & eos_name)
{
  tbqs(T_freeze_out_at_mu_eq_0, 0, 0, 0, eos_name);
  return eVal;
}

double EquationOfState::sfreeze(double T_freeze_out_at_mu_eq_0, const string & eos_name)
{
  return s_terms_T(T_freeze_out_at_mu_eq_0, eos_name);
}



////////////////////////////////////////////////////////////////////////////////
// FUNCTIONS TO UPDATE LOCATION IN PHASE DIAGRAM AND COMPUTE CORRESPONDING
// THERMODYNAMICS QUANTITIES NEEDED IN HYDRO


////////////////////////////////////////////////
// update phase diagram location given (s,B,S,Q)
bool EquationOfState::update_s(double sin) { return update_s(sin, 0.0, 0.0, 0.0); }
bool EquationOfState::update_s( double sin, double Bin, double Sin, double Qin,
                                bool verbose )
{
  print_now = verbose;

  bool success = false;
  if ( use_rootfinder )
    success = rootfinder_update_s(sin, Bin, Sin, Qin);
  else
  {
    std::cerr << __PRETTY_FUNCTION__ << ":" << __LINE__ << ": Option not supported!" << std::endl;
    exit(1);
  }

  return success;
}


////////////////////////////////////////////////
bool EquationOfState::find_root_with_seed(
      const string & e_or_s_mode, double e_or_s_in,
      double Bin, double Sin, double Qin,
      pEoS_base this_eos, vector<double> & result )
{
    const double hc = constants::hbarc_MeVfm;
    if ( VERBOSE > 2 || print_now )
    {
      std::cout << " --> currently trying " << this_eos->name
                << " EoS for solution..." << std::endl;
      std::cout << "     - point(" << e_or_s_mode << "): "
                << e_or_s_in << "   " << Bin << "   " << Sin << "   " << Qin
                << std::endl;
      std::cout << "     - seed: "
                << result[0]*hc << "   " << result[1]*hc << "   "
                << result[2]*hc << "   " << result[3]*hc << std::endl;
    }

    // if not energy, do entropy
    if ( e_or_s_mode != "energy" )
      return rootfinder.find_root(
              e_or_s_mode, e_or_s_in, Bin, Sin, Qin,
              this_eos->sBSQ, this_eos->tbqs_minima, this_eos->tbqs_maxima,
              result );
    else
      return rootfinder.find_root(
              e_or_s_mode, e_or_s_in, Bin, Sin, Qin,
              this_eos->eBSQ, this_eos->tbqs_minima, this_eos->tbqs_maxima,
              result );
}





////////////////////////////////////////////////
bool EquationOfState::rootfinder_update_s( double sin, double Bin,
                                           double Sin, double Qin )
{
  const double hc = constants::hbarc_MeVfm;

  // take sign of densities using lambda function
  auto sgn = [](double val) -> double { return (0.0 < val) - (val < 0.0); };

  bool solution_found = false;
  bool skipping_EoSs = true;  // start out skipping EoSs by default
  vector<double> result;

  string eos_currently_trying = "";

  //////////////////////////////////////////////////
  // try each EoS in turn
  for ( const auto & this_eos : chosen_EOSs )
  {
    ////////////////////////////////////////////////
    // skip EoSs where the particle has previously failed to find a solution
    if ( skip_failed_EoS && skipping_EoSs )
    {
      if ( this_eos->name != current_eos_name )
        continue;
      else
        skipping_EoSs = false;  //stop skipping at this point
    }

    eos_currently_trying = this_eos->name;

    ////////////////////////////////////////////////
    // try current location
    result = tbqsPosition;
    solution_found = find_root_with_seed( "entropy", sin, Bin, Sin, Qin,
                                          this_eos, result );


    ////////////////////////////////////////////////
    // if rootfinder fails, try a different seed
    if (!solution_found)
    {
      // try default seed at zero density
      result = vector<double>({tbqsPosition[0],0.0,0.0,0.0});
      solution_found = find_root_with_seed( "entropy", sin, Bin, Sin, Qin,
                                            this_eos, result );
    }


    ////////////////////////////////////////////////
    // if rootfinder fails, try a different seed
    if (!solution_found)
    {
      // try forced seed
      result = vector<double>({580.0/hc, sgn(Bin), sgn(Qin), sgn(Sin)});
//      result = vector<double>({std::min(2500.0/hc, (this_eos->tbqs_maxima)[0]),
//                                sgn(Bin), sgn(Qin), sgn(Sin)});
      solution_found = find_root_with_seed( "entropy", sin, Bin, Sin, Qin,
                                            this_eos, result );
    }


    ////////////////////////////////////////////////
    // if rootfinder fails, try another different seed
    if (!solution_found)
    {
      // try conformal diagonal seed
      auto conf_diag_EoS = std::dynamic_pointer_cast<EoS_conformal_diagonal>
                                  ( chosen_EOS_map["conformal_diagonal"] );

      result = conf_diag_EoS->get_tbqs_seed_from_sBSQ( sin, Bin, Sin, Qin );

      solution_found = find_root_with_seed( "entropy", sin, Bin, Sin, Qin,
                                            this_eos, result );
    }


    ////////////////////////////////////////////////
    // if all else fails and relevant flag is set,
    // ZERO CHARGE DENSITIES
    if (zero_unsolvable_charge_densities && !solution_found)
    {
      result = vector<double>({580.0/hc, 0.0, 0.0, 0.0});
//      result = vector<double>({std::min(2500.0/hc, (this_eos->tbqs_maxima)[0]),
//                                0.0, 0.0, 0.0});
      solution_found = find_root_with_seed( "entropy", sin, 0.0, 0.0, 0.0,
                                            this_eos, result );
    }


    ////////////////////////////////////////////////
    // Gap-region fallback: triggered ONLY when the table EoS could not
    // converge for this (s,rhoB,...) target.  Must run *inside* the loop so
    // it intercepts the failure before the loop falls through to the
    // conformal EoSs (which are globally invertible and would otherwise
    // mask the gap and make this fallback dead code).
    if (!solution_found && this_eos->name == "table")
    {
      auto table_eos = std::dynamic_pointer_cast<EoS_table>(this_eos);
      if (table_eos && (table_eos->gap_table.loaded
                        || settingsPtr->gap_analytic_enabled))
      {
        double T_g = 0.0, muB_g = 0.0, muQ_g = 0.0, muS_g = 0.0;
        double thermo_g[N_GAP_THERMO];
        bool gap_hit = false;
        if (settingsPtr->gap_analytic_enabled)
          gap_hit = table_eos->gap_table.analytic_lookup(sin, Bin, Sin, Qin,
                      T_g, muB_g, muQ_g, muS_g, thermo_g);
        if (!gap_hit && table_eos->gap_table.loaded)
          gap_hit = table_eos->gap_table.lookup(sin, Bin, Sin, Qin,
                      settingsPtr->gap_lookup_mode,
                      settingsPtr->gap_match_tolerance,
                      T_g, muB_g, muQ_g, muS_g, thermo_g, /*use_energy=*/false);
        if (gap_hit)
        {
          tbqs_from_gap(T_g, muB_g, muQ_g, muS_g, thermo_g);
          solution_found = true;
          if (VERBOSE > 0)
            std::cout << " --> gap fallback (entropy): (T,muB)="
                      << T_g << "," << muB_g << std::endl;
          break;  // gap hit IS the table-EoS solution; do not advance to conformal
        }
      }
    }


    ////////////////////////////////////////////////
    // stop iterating through available EoSs when solution found
    if (solution_found)
    {
      tbqs( result, this_eos ); // set thermodynamics using solution

      if ( VERBOSE > 2 || print_now )
      {
        std::cout << " --> found a solution with "
                  << this_eos->name << " EoS!" << std::endl;
        std::cout << " --> solution has (T,muB,muS,muQ) = "
                  << T() << "   " << muB() << "   " << muS() << "   " << muQ()
                  << std::endl;
        std::cout << " --> solution has (s,B,S,Q) = "
                  << s() << "   " << B() << "   " << S() << "   " << Q()
                  << std::endl;
      }

      //========================================================================
      // check if cs2 or the mu/T ratios are going haywire, in which case,
      // don't trust this EoS!  Move on to the next one instead...
      if ( settingsPtr->prohibit_unstable_cs2 && cs2() < 0.0 )
        continue;
      else if ( settingsPtr->prohibit_acausal_cs2 && cs2() > 1.0 )
        continue;
//      else if ( restrict_mu_T_ratios && this_eos->name == "table"
//                && sqrt(muB()*muB()+muS()*muS()+muQ()*muQ()) > 4.0*T() )
      else if ( settingsPtr->restrict_mu_T_ratios && this_eos->name == "table"
                && std::max( std::max( std::abs(muB()), std::abs(muS()) ),
                        std::abs(muQ()) ) > 3.5*T() )
        continue;
      else
        break;
    }
  }

  if (!solution_found)
  {
    std::cout << "No solution found!" << std::endl;
    std::cout << "Last attempted EoS: " << eos_currently_trying << endl;
    std::cerr << "No solution found!" << std::endl;
    std::cout << "Last attempted EoS: " << eos_currently_trying << endl;
    std::cout << "Failed to find a solution for (s,B,S,Q) = "
              << sin << "   " << Bin << "   " << Sin << "   " << Qin << std::endl;
    std::stringstream ss;
    ss << "(s,B,S,Q) = " << sin << "   " << Bin << "   " << Sin << "   " << Qin;
    throw std::invalid_argument(ss.str());
  }

  return solution_found;
}
////////////////////////////////////////////////



////////////////////////////////////////////////
// update phase diagram location given (e,B,S,Q) and return resulting s
double EquationOfState::s_out( double ein, bool & solution_found )
                        { return s_out(ein, 0.0, 0.0, 0.0, solution_found); }
double EquationOfState::s_out( double ein, double Bin, double Sin,
                               double Qin, bool & solution_found, bool verbose )
{
  print_now = verbose;

  double result = 0.0;
  if ( use_rootfinder )
    result = rootfinder_s_out(ein, Bin, Sin, Qin, solution_found);
  else
  {
    std::cerr << __PRETTY_FUNCTION__ << ":" << __LINE__ << ": Option not supported!" << std::endl;
    exit(1);
  }

  return (result);
}

////////////////////////////////////////////////////////////////////////////////
double EquationOfState::rootfinder_s_out( double ein, double Bin, double Sin,
                                          double Qin, bool & solution_found )
{
  const double hc = constants::hbarc_MeVfm;
  std::string eos_type;
  // take sign of densities using lambda function
  auto sgn = [](double val) -> double { return (0.0 < val) - (val < 0.0); };

  // used for seed value in rootfinder
  vector<double> result;
  solution_found = false;

  ///////////////////////////////////////////////////////////
  // try each EoS in turn
  for ( const auto & this_eos : chosen_EOSs )
  {
    //assign eos type string
    eos_type = this_eos->name;
    /////////////////////////////////////////////////////////
    // try forced seed first
    result = vector<double>({580.0/hc, sgn(Bin), sgn(Qin), sgn(Sin)});
//    result = vector<double>({std::min(2500.0/hc, (this_eos->tbqs_maxima)[0]),
//                                sgn(Bin), sgn(Qin), sgn(Sin)});
    solution_found = find_root_with_seed( "energy", ein, Bin, Sin, Qin,
                                          this_eos, result );


    /////////////////////////////////////////////////////////
    // try conformal diagonal seed next
    if (!solution_found)
    {
      // try conformal diagonal seed
      auto conf_diag_EoS = std::dynamic_pointer_cast<EoS_conformal_diagonal>
                                  ( chosen_EOS_map["conformal_diagonal"] );

      result = conf_diag_EoS->get_tbqs_seed_from_eBSQ( ein, Bin, Sin, Qin );

      solution_found = find_root_with_seed( "energy", ein, Bin, Sin, Qin,
                                            this_eos, result );
    }

    ////////////////////////////////////////////////
    // if all else fails and relevant flag is set,
    // ZERO CHARGE DENSITIES
    if (zero_unsolvable_charge_densities && !solution_found)
    {
      result = vector<double>({580.0/hc, 0.0, 0.0, 0.0});
//      result = vector<double>({std::min(2500.0/hc, (this_eos->tbqs_maxima)[0]),
//                                0.0, 0.0, 0.0});
      solution_found = find_root_with_seed( "energy", ein, 0.0, 0.0, 0.0,
                                            this_eos, result );
    }


    /////////////////////////////////////////////////////////
    // Gap-region fallback (energy path): triggered ONLY when the table EoS
    // failed to converge.  Must run inside the loop so it intercepts the
    // failure before falling through to the conformal EoSs.
    if (!solution_found && this_eos->name == "table")
    {
      auto table_eos = std::dynamic_pointer_cast<EoS_table>(this_eos);
      if (table_eos && (table_eos->gap_table.loaded
                        || settingsPtr->gap_analytic_enabled))
      {
        double T_g = 0.0, muB_g = 0.0, muQ_g = 0.0, muS_g = 0.0;
        double thermo_g[N_GAP_THERMO];
        bool gap_hit = false;
        if (settingsPtr->gap_analytic_enabled)
          gap_hit = table_eos->gap_table.analytic_lookup(ein, Bin, Sin, Qin,
                      T_g, muB_g, muQ_g, muS_g, thermo_g);
        if (!gap_hit && table_eos->gap_table.loaded)
          gap_hit = table_eos->gap_table.lookup(ein, Bin, Sin, Qin,
                      settingsPtr->gap_lookup_mode,
                      settingsPtr->gap_match_tolerance,
                      T_g, muB_g, muQ_g, muS_g, thermo_g, /*use_energy=*/true);
        if (gap_hit)
        {
          tbqs_from_gap(T_g, muB_g, muQ_g, muS_g, thermo_g);
          solution_found = true;
          if (VERBOSE > 0)
            std::cout << " --> gap fallback (energy): (T,muB)="
                      << T_g << "," << muB_g << std::endl;
          break;  // gap hit IS the table-EoS solution; do not advance to conformal
        }
      }
    }


    /////////////////////////////////////////////////////////
    // stop iterating through available EoSs when solution found
    if (solution_found)
    {
      tbqs( result, this_eos ); // set thermodynamics using solution

      if ( VERBOSE > 2 || print_now )
      {
        std::cout << " --> found a solution with "
                  << this_eos->name << " EoS!" << std::endl;
        std::cout << " --> solution has (T,muB,muS,muQ) = "
                  << T() << "   " << muB() << "   " << muS() << "   " << muQ()
                  << std::endl;
        std::cout << " --> solution has (e,B,S,Q) = "
                  << e() << "   " << B() << "   " << S() << "   " << Q()
                  << std::endl;
      }

      //========================================================================
      // check if cs2 or the mu/T ratios are going haywire, in which case,
      // don't trust this EoS!  Move on to the next one instead...
      if ( settingsPtr->prohibit_unstable_cs2 && cs2() < 0.0 )
        continue;
      else if ( settingsPtr->prohibit_acausal_cs2 && cs2() > 1.0 )
        continue;
//      else if ( restrict_mu_T_ratios && this_eos->name == "table"
//                && sqrt(muB()*muB()+muS()*muS()+muQ()*muQ()) > 4.0*T() )
      else if ( settingsPtr->restrict_mu_T_ratios && this_eos->name == "table"
                && std::max( std::max( std::abs(muB()), std::abs(muS()) ),
                        std::abs(muQ()) ) > 3.5*T() )
        continue;
      else
        break;
    }
  }

  if (!solution_found)
  {
    //check this eos type
    std::cout << "No solution found!" << std::endl;
    //std::cerr << "No solution found!" << std::endl;
    bool test = eBSQ_has_solution_in_conformal_diagonal(ein, Bin, Sin, Qin);
    std::cout << "Conformal diagonal EoS should have sol: " << test << std::endl;
    std::cout << "Last attempted EoS: " << eos_type << endl;
    std::cout << "Failed to find a solution for (e,B,S,Q) = "
              << ein << "   " << Bin << "   " << Sin << "   " << Qin << std::endl;
    std::stringstream ss;
    ss << "(e,B,S,Q) = " << ein << "   " << Bin << "   " << Sin << "   " << Qin;
    //throw std::invalid_argument(ss.str());
  }

//  cout << "Exiting " << __FUNCTION__ << endl;

  // this is set in most recent call to tbqs()
  return entrVal;

}
////////////////////////////////////////////////////////////////////////////////

