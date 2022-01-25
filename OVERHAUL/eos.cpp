#include "eos.h"

InterpolatorND<4> EquationOfState::equation_of_state_table;

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

#include "eos_derivatives.cpp"
#include "eos_initialization.cpp"

#include "constants.h"

using namespace constants;

using std::vector;
using std::string;

////////////////////////////////////////////////////////////////////////////////
// Constructors
EquationOfState::EquationOfState(){}

EquationOfState::EquationOfState(string quantityFile, string derivFile)
{
    init(quantityFile, derivFile);
}



////////////////////////////////////////////////////////////////////////////////
bool EquationOfState::point_not_in_range(
                        double setT, double setmuB, double setmuQ,
                        double setmuS, bool use_conformal )
{
  auto & tbqs_mins = ( use_conformal ) ? conformal_tbqs_minima : tbqs_minima;
  auto & tbqs_maxs = ( use_conformal ) ? conformal_tbqs_maxima : tbqs_maxima;
  string conformal_string = ( use_conformal ) ? "conformal" : "non-conformal";

  if(setT < tbqs_mins[0] || setT > tbqs_maxs[0])
  { 
    std::cout << "T = " << setT
      << " is out of (" << conformal_string << ") range."
         " Valid values are between ["
      << tbqs_mins[0] << "," << tbqs_maxs[0] << "]" << std::endl;
    return true;
  }
  if(setmuB < tbqs_mins[1] || setmuB > tbqs_maxs[1])
  {
    std::cout << "muB = " << setmuB
      << " is out of (" << conformal_string << ") range."
         " Valid values are between ["
      << tbqs_mins[1] << "," << tbqs_maxs[1] << "]" << std::endl;
    return true;
  }
  if(setmuQ < tbqs_mins[2] || setmuQ > tbqs_maxs[2])
  {
    std::cout << "muQ = " << setmuQ
      << " is out of (" << conformal_string << ") range."
         " Valid values are between ["
      << tbqs_mins[2] << "," << tbqs_maxs[2] << "]" << std::endl;
    return true;
  }
  if(setmuS < tbqs_mins[3] || setmuS > tbqs_maxs[3])
  {
    std::cout << "muS = " << setmuS
      << " is out of (" << conformal_string << ") range."
         " Valid values are between ["
      << tbqs_mins[3] << "," << tbqs_maxs[3] << "]" << std::endl;
    return true;
  }
  return false;
}


////////////////////////////////////////////////////////////////////////////////
void EquationOfState::tbqs( vector<double> & tbqsIn, bool use_conformal )
{
  tbqs( tbqsIn[0], tbqsIn[1], tbqsIn[2], tbqsIn[3], use_conformal);
}

////////////////////////////////////////////////////////////////////////////////
void EquationOfState::tbqs( double setT, double setmuB, double setmuQ,
                            double setmuS, bool use_conformal )
{
  if ( point_not_in_range( setT, setmuB, setmuQ, setmuS, use_conformal ) )
    return;

  tbqsPosition[0] = setT;
	tbqsPosition[1] = setmuB;
	tbqsPosition[2] = setmuQ;
	tbqsPosition[3] = setmuS;

  // if we are in range, compute all thermodynamic quantities at the new point
  evaluate_thermodynamics(use_conformal);
}


void EquationOfState::evaluate_thermodynamics(bool use_conformal)
{
  if ( settingsPtr->EoS_type == "Conformal" or use_conformal )
  {
    // EXPECTS UNITS OF MEV!!!
    double phase_diagram_point[4]
        = { tbqsPosition[0], tbqsPosition[1], tbqsPosition[2], tbqsPosition[3] };

    double thermodynamics[17];
    eos_conformal::get_full_thermo(phase_diagram_point, thermodynamics);

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

    /*for (int i = 0; i < 4; i++)
      cout << tbqsPosition[i] << "   ";
    for (int i = 0; i < 17; i++)
      cout << thermodynamics[i] << "   ";
    cout << endl;*/
  }
  else if ( use_static_C_library )
  {
    // EXPECTS UNITS OF MEV!!!
    double phase_diagram_point[4]
        = { tbqsPosition[0], tbqsPosition[1], tbqsPosition[2], tbqsPosition[3] };

/*phase_diagram_point[0] = 150.0;
phase_diagram_point[1] = 0.0;
phase_diagram_point[2] = 0.0;
phase_diagram_point[3] = 0.0;*/

    double thermodynamics[17];
    STANDARD_get_full_thermo(phase_diagram_point, thermodynamics);

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
  }
  else
  {
    vector<double> thermodynamics;  // gets re-sized inside evaluate function

    // evaluate EoS interpolator at current location (S and Q NOT SWAPPED)
    equation_of_state_table.evaluate( tbqsPosition, thermodynamics ); 

    if ( thermodynamics.size() != 17 )
    {
      cerr << "PROBLEM" << endl;
      exit(1);
    }

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
  }

/*cout << "THERMO DUMP: " << pVal << "   " << entrVal << "   " << BVal << "   "
      << SVal << "   " << QVal << "   " << eVal << "   " << cs2Val << "   "
      << db2 << "   " << dq2 << "   " << ds2 << "   " << dbdq << "   "
      << dbds << "   " << dsdq << "   " << dtdb << "   " << dtdq << "   "
      << dtds << "   " << dt2 << endl;

if (true) exit(1);*/
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


double EquationOfState::dwds()
{
	double charge_terms	/*if charge densities are not all zero*/
			= ( abs(BVal)>1e-10 || abs(SVal)>1e-10 || abs(QVal)>1e-10 ) ?
			  BVal/dentr_dmub() + QVal/dentr_dmuq() + SVal/dentr_dmus() : 0.0;

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
	double charge_terms	/*if charge densities are not all zero*/
			= ( abs(BVal)>1e-10 || abs(SVal)>1e-10 || abs(QVal)>1e-10 ) ?
			  entrVal/db_dt() + BVal/db_dmub() + QVal/db_dmuq() + SVal/db_dmus() : 0.0;

    return muB() + charge_terms;
}

double EquationOfState::dwdS()
{
	double charge_terms	/*if charge densities are not all zero*/
			= ( abs(BVal)>1e-10 || abs(SVal)>1e-10 || abs(QVal)>1e-10 ) ?
			  entrVal/ds_dt() + BVal/ds_dmub() + QVal/ds_dmuq() + SVal/ds_dmus() : 0.0;

    return muS() + charge_terms;
}

double EquationOfState::dwdQ()
{
	double charge_terms	/*if charge densities are not all zero*/
			= ( abs(BVal)>1e-10 || abs(SVal)>1e-10 || abs(QVal)>1e-10 ) ?
			  entrVal/dq_dt() + BVal/dq_dmub() + QVal/dq_dmuq() + SVal/dq_dmus() : 0.0;

    return muQ() + charge_terms;
}


double EquationOfState::cs2out(double Tt) {  //return cs2 given t and mu's=0
cout << "Entering " << __PRETTY_FUNCTION__ << "::" << __LINE__ << endl;
    tbqs(Tt, 0.0, 0.0, 0.0, false);
cout << "Leaving " << __PRETTY_FUNCTION__ << "::" << __LINE__ << endl;
    return cs2Val;
}

double EquationOfState::cs2out(double Tt, double muBin, double muQin, double muSin) {  //return cs2 given t and mu's
cout << "Entering " << __PRETTY_FUNCTION__ << "::" << __LINE__ << endl;
    tbqs(Tt, muBin, muQin, muSin, false);
cout << "Leaving " << __PRETTY_FUNCTION__ << "::" << __LINE__ << endl;
    return cs2Val;
}

double EquationOfState::wfz(double Tt) {   // return e + p for tbqs
cout << "Entering " << __PRETTY_FUNCTION__ << "::" << __LINE__ << endl;
    tbqs(Tt, 0.0, 0.0, 0.0, false);
cout << "Leaving " << __PRETTY_FUNCTION__ << "::" << __LINE__ << endl;
    return eVal + pVal;
}

double EquationOfState::wfz(double Tt, double muBin, double muQin, double muSin) {   // return e + p for tbqs
cout << "Entering " << __PRETTY_FUNCTION__ << "::" << __LINE__ << endl;
    tbqs(Tt, muBin, muQin, muSin, false);
cout << "Leaving " << __PRETTY_FUNCTION__ << "::" << __LINE__ << endl;
    return eVal + pVal;
}


double EquationOfState::s_terms_T(double Tt)
{
cout << "Entering " << __PRETTY_FUNCTION__ << "::" << __LINE__ << endl;
  tbqs(Tt, 0, 0, 0, false);
cout << "Leaving " << __PRETTY_FUNCTION__ << "::" << __LINE__ << endl;
  return entrVal;
}


// UNCOMMENTED BY C. PLUMBERG
void EquationOfState::eosin(std::string type){}

double EquationOfState::A() { return w()-s()*dwds(); }


// confirm with Jaki
double EquationOfState::efreeze(double T_freeze_out_at_mu_eq_0)
{
  tbqs(T_freeze_out_at_mu_eq_0, 0, 0, 0, false);
  return eVal;
}

double EquationOfState::sfreeze(double T_freeze_out_at_mu_eq_0)
{
  return s_terms_T(T_freeze_out_at_mu_eq_0, false);
}



////////////////////////////////////////////////////////////////////////////////
// FUNCTIONS TO UPDATE LOCATION IN PHASE DIAGRAM AND COMPUTE CORRESPONDING
// THERMODYNAMICS QUANTITIES NEEDED IN HYDRO
// USE EITHER ROOTFINDER OR DELAUNAY INTERPOLATION


////////////////////////////////////////////////
// update phase diagram location given (s,B,S,Q)
bool EquationOfState::update_s(double sin) { return update_s(sin, 0.0, 0.0, 0.0); }
bool EquationOfState::update_s(double sin, double Bin, double Sin, double Qin)
{
  bool success = false;
  if ( use_delaunay )
    success = delaunay_update_s(sin, Bin, Sin, Qin);
  else if ( use_rootfinder )
    success = rootfinder_update_s(sin, Bin, Sin, Qin);
  else
  {
    std::cerr << __PRETTY_FUNCTION__ << ":" << __LINE__ << ": Option not supported!" << std::endl;
    exit(1);
  }

  return success;
}

////////////////////////////////////////////////
bool EquationOfState::delaunay_update_s(double sin, double Bin, double Sin, double Qin)
{
  if (true)
  {
    std::cout << "You still need to check units!  Etc." << std::endl;
    std::cerr << "You still need to check units!  Etc." << std::endl;
    exit(1);
  }
  vector<double> result = tbqsPosition;

  //bool success = entr_delaunay.interpolate( {sin, Bin, Sin, Qin}, result, true );
  bool success = false;
  tbqs( result, false );
  return success;
}

////////////////////////////////////////////////
bool EquationOfState::rootfinder_update_s(double sin, double Bin, double Sin, double Qin)
{
  vector<double> result = tbqsPosition;
  
  bool success = rootfinder.find_sBSQ_root( sin, Bin, Sin, Qin, sBSQ_functional,
                                            tbqs_minima, tbqs_maxima, result );
  if (success)
  {
    tbqs( result, false );
    return success;
  }
  else
  {
    success = rootfinder.find_sBSQ_root( sin, Bin, Sin, Qin, conformal_sBSQ_functional,
                                         conformal_tbqs_minima, conformal_tbqs_maxima, result );
    tbqs( result, true );
    return success;
  }
}
////////////////////////////////////////////////



////////////////////////////////////////////////
// update phase diagram location given (e,B,S,Q) and return resulting s
double EquationOfState::s_out(double ein) { return s_out(ein, 0.0, 0.0, 0.0); }
double EquationOfState::s_out(double ein, double Bin, double Sin, double Qin)
{
  double result = 0.0;
  if ( use_delaunay )
    result = delaunay_s_out(ein, Bin, Sin, Qin);
  else if ( use_rootfinder )
    result = rootfinder_s_out(ein, Bin, Sin, Qin);
  else
  {
    std::cerr << __PRETTY_FUNCTION__ << ":" << __LINE__ << ": Option not supported!" << std::endl;
    exit(1);
  }

  return (result);
}

////////////////////////////////////////////////
double EquationOfState::delaunay_s_out(double ein, double Bin, double Sin, double Qin)
{
  if (true)
  {
    std::cout << "You still need to check units!  Etc." << std::endl;
    std::cerr << "You still need to check units!  Etc." << std::endl;
    exit(1);
  }
  vector<double> result = tbqsPosition;
  //e_delaunay.interpolate( {ein, Bin, Sin, Qin}, result, true );
  tbqs( result, false );
  return entrVal;
}

////////////////////////////////////////////////
double EquationOfState::rootfinder_s_out(double ein, double Bin, double Sin, double Qin)
{
  vector<double> result = tbqsPosition;
  bool success = rootfinder.find_eBSQ_root( ein, Bin, Sin, Qin, eBSQ_functional,
                                            tbqs_minima, tbqs_maxima, result );
  if (success)
  {
    tbqs( result, false );
    return entrVal;
  }
  else
  {
    success = rootfinder.find_eBSQ_root( ein, Bin, Sin, Qin, conformal_eBSQ_functional,
                                         conformal_tbqs_minima, conformal_tbqs_maxima, result );
    tbqs( result, true );
    return entrVal;
  }
}
////////////////////////////////////////////////


// point     = (T,muB,muQ,muS)
// densities = (e,rhoB,rhoS,rhoQ)
void EquationOfState::get_eBSQ_densities_from_interpolator(
        double point[], double densities[] )  // point and densities both length = 4
{
    vector<double> results;
    equation_of_state_table.evaluate(
      vector<double>(point, point + 4), results,
      vector<string>({ "e","B","S","Q" }) );
    std::copy(results.begin(), results.end(), densities);
}

// point     = (T,muB,muQ,muS)
// densities = (s,rhoB,rhoS,rhoQ)
void EquationOfState::get_sBSQ_densities_from_interpolator(
        double point[], double densities[] )  // point and densities both length = 4
{
    vector<double> results;
    equation_of_state_table.evaluate(
      vector<double>(point, point + 4), results,
      vector<string>({ "s","B","S","Q" }) );
    std::copy(results.begin(), results.end(), densities);
}