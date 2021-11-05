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
#include "../eos_delaunay/eos_delaunay.h"

#include "../constants.h"

using namespace constants;

using std::vector;
using std::string;

constexpr bool use_exact = true;
constexpr bool accept_nearest_neighbor = false;
constexpr bool discard_unsolvable_charge_densities = false;
//constexpr bool check_derivatives = false;

constexpr size_t STEPS = 1000000;
constexpr int VERBOSE = 0;
constexpr double TOLERANCE = 1e-12;

////////////////////////////////////////////////////////////////////////////////
// Constructors
EquationOfState::EquationOfState(){}

EquationOfState::EquationOfState(string quantityFile, string derivFile)
{
    init(quantityFile, derivFile);
}



////////////////////////////////////////////////////////////////////////////////
void EquationOfState::tbqs( vector<double> & tbqsIn )
{
  tbqs( tbqsIn[0], tbqsIn[1], tbqsIn[2], tbqsIn[3] );
}

////////////////////////////////////////////////////////////////////////////////
void EquationOfState::tbqs(double setT, double setmuB, double setmuQ, double setmuS)
{
  if(setT < minT || setT > maxT)
  {
    std::cout << "T = " << setT << " is out of range. Valid values are between ["
      << minT << "," << maxT << "]" << std::endl;
    return;
  }
  if(setmuB < minMuB || setmuB > maxMuB)
  {
    std::cout << "muB = " << setmuB << " is out of range. Valid values are between ["
      << minMuB << "," << maxMuB << "]" << std::endl;
    return;
  }
  if(setmuQ < minMuQ || setmuQ > maxMuQ)
  {
    std::cout << "muQ = " << setmuQ << " is out of range. Valid values are between ["
      << minMuQ << "," << maxMuQ << "]" << std::endl;
    return;
  }
  if(setmuS < minMuS || setmuS > maxMuS)
  {
    std::cout << "muS = " << setmuS << " is out of range. Valid values are between ["
      << minMuS << "," << maxMuS << "]" << std::endl;
    return;
  }

  tbqsPosition[0] = setT;
	tbqsPosition[1] = setmuB;
	tbqsPosition[2] = setmuQ;
	tbqsPosition[3] = setmuS;

  // if we are in range, compute all thermodynamic quantities at the new point
  evaluate_thermodynamics();
}


void EquationOfState::evaluate_thermodynamics()
{
	// EXPECTS UNITS OF MEV!!!
	double phase_diagram_point[4]	// NOTE: S <<-->> Q swapped!!!
			= {setT*hbarc_MeVfm, setmuB*hbarc_MeVfm, setmuS*hbarc_MeVfm, setmuQ*hbarc_MeVfm};
	double thermodynamics[17];
  get_full_thermo(phase_diagram_point, thermodynamics);

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
    tbqs(Tt, 0.0, 0.0, 0.0);
    return cs2Val;
}

double EquationOfState::cs2out(double Tt, double muBin, double muQin, double muSin) {  //return cs2 given t and mu's
    tbqs(Tt, muBin, muQin, muSin);
    return cs2Val;
}

double EquationOfState::wfz(double Tt) {   // return e + p for tbqs
    tbqs(Tt, 0.0, 0.0, 0.0);
    return eVal + pVal;
}

double EquationOfState::wfz(double Tt, double muBin, double muQin, double muSin) {   // return e + p for tbqs
    tbqs(Tt, muBin, muQin, muSin);
    return eVal + pVal;
}


double EquationOfState::s_terms_T(double Tt)
{
  tbqs(Tt, 0, 0, 0);
  return entrVal;
}


// UNCOMMENTED BY C. PLUMBERG
void EquationOfState::eosin(std::string type){}

double EquationOfState::A() { return w()-s()*dwds(); }


// confirm with Jaki
double EquationOfState::efreeze(double T_freeze_out_at_mu_eq_0)
{
  tbqs(T_freeze_out_at_mu_eq_0, 0, 0, 0);
  return eVal;
}

double EquationOfState::sfreeze(double T_freeze_out_at_mu_eq_0)
{
  return s_terms_T(T_freeze_out_at_mu_eq_0);
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
  vector<double> result(4, 0.0);

  //bool success = entr_delaunay.interpolate( {sin, Bin, Sin, Qin}, result, true );
  bool success = false;
  tbqs( result );
  return success;
}

////////////////////////////////////////////////
bool EquationOfState::rootfinder_update_s(double sin, double Bin, double Sin, double Qin)
{
  vector<double> result(4, 0.0);
  
  bool success = rootfinder.find_sBSQ_root( sin, Bin, Sin, Qin, sBSQ_functional, result );
  tbqs( result );
  return success;
}
////////////////////////////////////////////////



////////////////////////////////////////////////
// update phase diagram location given (e,B,S,Q) and return resulting s
double EquationOfState::s_out(double ein) { return s_out(ein, 0.0, 0.0, 0.0); }
double EquationOfState::s_out(double ein, double Bin, double Sin, double Qin)
{
  double result = 0.0;
  if ( use_delaunay )
    result = delaunay_s_out(sin, Bin, Sin, Qin);
  else if ( use_rootfinder )
    result = rootfinder_s_out(sin, Bin, Sin, Qin);
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
  vector<double> result(4, 0.0);
  //e_delaunay.interpolate( {ein, Bin, Sin, Qin}, result, true );
  tbqs( result );
  return entrVal;
}

////////////////////////////////////////////////
double EquationOfState::rootfinder_s_out(double ein, double Bin, double Sin, double Qin)
{
  vector<double> result(4, 0.0);
  rootfinder.find_eBSQ_root( ein, Bin, Sin, Qin, eBSQ_functional, result );
  tbqs( result );
  return entrVal;
}
////////////////////////////////////////////////