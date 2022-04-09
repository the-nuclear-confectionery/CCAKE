//#include "eos.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "constants.h"
#include "eos.h"
#include "eos_delaunay/eos_delaunay.h"

#include "eos_derivatives.cpp"
#include "eos_initialization.cpp"

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
      { if ( eos_name == "default" )
          tbqs( setT, setmuB, setmuQ, setmuS, chosen_EOS_map[default_eos_name] );
        else
          tbqs( setT, setmuB, setmuQ, setmuS, chosen_EOS_map[eos_name] ); }

void EquationOfState::tbqs( vector<double> & tbqsIn, const string & eos_name )
      { if ( eos_name == "default" )
          tbqs( tbqsIn, chosen_EOS_map[default_eos_name] );
        else
          tbqs( tbqsIn, chosen_EOS_map[eos_name] ); }



////////////////////////////////////////////////////////////////////////////////
void EquationOfState::tbqs( double setT, double setmuB, double setmuQ,
                            double setmuS, pEoS_base peos )
{
//cout << __PRETTY_FUNCTION__ << ": " << peos->name << endl;
  bool point_is_in_range = !point_not_in_range( setT, setmuB, setmuQ, setmuS, peos );
  if ( point_is_in_range )
  {
    tbqsPosition[0] = setT;
    tbqsPosition[1] = setmuB;
    tbqsPosition[2] = setmuQ;
    tbqsPosition[3] = setmuS;

    // if we are in range, compute all thermodynamic quantities at the new point
    evaluate_thermodynamics( peos );
  }
  else
  {
    std::cout << "Point was out of range! You need to tell me what to do!" << std::endl;
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


double EquationOfState::cs2out(double Tt, const string & eos_name)
{  //return cs2 given t and mu's=0
  std::cout << __PRETTY_FUNCTION__ << "::" << __LINE__ << ": "
            << "using " << eos_name << "!" << std::endl;
  tbqs(Tt, 0.0, 0.0, 0.0, eos_name);
  return cs2Val;
}

double EquationOfState::cs2out(double Tt, double muBin, double muQin, double muSin, const string & eos_name)
{  //return cs2 given t and mu's
  std::cout << __PRETTY_FUNCTION__ << "::" << __LINE__ << ": "
            << "using " << eos_name << "!" << std::endl;
  tbqs(Tt, muBin, muQin, muSin, eos_name);
  return cs2Val;
}

double EquationOfState::wfz(double Tt, const string & eos_name)
{   // return e + p for tbqs
  std::cout << __PRETTY_FUNCTION__ << "::" << __LINE__ << ": "
            << "using " << eos_name << "!" << std::endl;
  tbqs(Tt, 0.0, 0.0, 0.0, eos_name);
  return eVal + pVal;
}

double EquationOfState::wfz(double Tt, double muBin, double muQin, double muSin, const string & eos_name)
{   // return e + p for tbqs
  std::cout << __PRETTY_FUNCTION__ << "::" << __LINE__ << ": "
            << "using " << eos_name << "!" << std::endl;
  tbqs(Tt, muBin, muQin, muSin, eos_name);
  return eVal + pVal;
}


double EquationOfState::s_terms_T(double Tt, const string & eos_name)
{
  std::cout << __PRETTY_FUNCTION__ << "::" << __LINE__ << ": "
            << "using " << eos_name << "!" << std::endl;
  tbqs(Tt, 0, 0, 0, eos_name);
  return entrVal;
}


// UNCOMMENTED BY C. PLUMBERG
void EquationOfState::eosin(std::string type){}

double EquationOfState::A() { return w()-s()*dwds(); }


// confirm with Jaki
double EquationOfState::efreeze(double T_freeze_out_at_mu_eq_0, const string & eos_name)
{
  std::cout << __PRETTY_FUNCTION__ << "::" << __LINE__ << ": "
            << "using " << eos_name << "!" << std::endl;
  tbqs(T_freeze_out_at_mu_eq_0, 0, 0, 0, eos_name);
cout << __PRETTY_FUNCTION__ << "::" << __LINE__ << endl;
  return eVal;
}

double EquationOfState::sfreeze(double T_freeze_out_at_mu_eq_0, const string & eos_name)
{
  return s_terms_T(T_freeze_out_at_mu_eq_0, eos_name);
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
  tbqs( result, "table" );
  return success;
}

////////////////////////////////////////////////
bool EquationOfState::rootfinder_update_s(double sin, double Bin, double Sin, double Qin)
{
  bool success = false;
  vector<double> result = tbqsPosition;

  // try each EoS in turn
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  for ( const auto & this_eos : chosen_EOSs )
  {
    std::cout << " --> currently trying " << this_eos->name << " EoS for solution..." << std::endl;
    result = tbqsPosition;
    std::cout << "     - seed: "
              << result[0] << "   " << result[1] << "   "
              << result[2] << "   " << result[3] << std::endl;
    success
      = rootfinder.find_sBSQ_root( sin, Bin, Sin, Qin, this_eos->sBSQ,
                                   this_eos->tbqs_minima, this_eos->tbqs_maxima,
                                   result );

    // try a different seed value if default guess fails
    if (!success)
    {
      // try twice the grid maxima
      result = this_eos->get_tbqs_maxima_no_ext();
      std::for_each(result.begin(), result.end(), [](double &c){ c *= 2.0; });
      std::cout << "     - seed: "
                << result[0] << "   " << result[1] << "   "
                << result[2] << "   " << result[3] << std::endl;
      success
        = rootfinder.find_sBSQ_root( sin, Bin, Sin, Qin, this_eos->sBSQ,
                                     this_eos->tbqs_minima, this_eos->tbqs_maxima,
                                     result );
    }

    // stop iterating through available EoSs when solution found
    if (success)
    {
      std::cout << " --> found a solution with " << this_eos->name << " EoS!" << std::endl;
      current_eos_name = this_eos->name;
      tbqs( result, this_eos ); // set thermodynamics using solution
      break;
    }
  }

  if (!success)
  {
    std::cout << "No solution found!" << std::endl;
    std::cerr << "No solution found!" << std::endl;
    exit(101);
  }

  return success;
}
////////////////////////////////////////////////



////////////////////////////////////////////////
// update phase diagram location given (e,B,S,Q) and return resulting s
double EquationOfState::s_out( double ein, bool & solution_found )
                        { return s_out(ein, 0.0, 0.0, 0.0, solution_found); }
double EquationOfState::s_out( double ein, double Bin, double Sin,
                               double Qin, bool & solution_found )
{
  double result = 0.0;
  if ( use_delaunay )
    result = delaunay_s_out(ein, Bin, Sin, Qin, solution_found);
  else if ( use_rootfinder )
    result = rootfinder_s_out(ein, Bin, Sin, Qin, solution_found);
  else
  {
    std::cerr << __PRETTY_FUNCTION__ << ":" << __LINE__ << ": Option not supported!" << std::endl;
    exit(1);
  }

  return (result);
}

////////////////////////////////////////////////
double EquationOfState::delaunay_s_out( double ein, double Bin, double Sin,
                                        double Qin, bool & solution_found )
{
  if (true)
  {
    std::cout << "You still need to check units!  Etc." << std::endl;
    std::cerr << "You still need to check units!  Etc." << std::endl;
    exit(1);
  }
  vector<double> result = tbqsPosition;
  //e_delaunay.interpolate( {ein, Bin, Sin, Qin}, result, true );
  tbqs( result, "table" );
  return entrVal;
}

////////////////////////////////////////////////////////////////////////////////
double EquationOfState::rootfinder_s_out( double ein, double Bin, double Sin,
                                          double Qin, bool & solution_found )
{
  // used for seed value in rootfinder
  vector<double> result = tbqsPosition;

  // try each EoS in turn
  std::cout << __PRETTY_FUNCTION__ << std::endl;
  for ( const auto & this_eos : chosen_EOSs )
  {
    std::cout << " --> currently trying " << this_eos->name << " EoS for solution..." << std::endl;
    result = tbqsPosition;
    std::cout << "     - seed: "
              << result[0]*hc << "   " << result[1]*hc << "   "
              << result[2]*hc << "   " << result[3]*hc << std::endl;
    solution_found
      = rootfinder.find_eBSQ_root( ein, Bin, Sin, Qin, this_eos->eBSQ,
                                   this_eos->tbqs_minima, this_eos->tbqs_maxima,
                                   result );

    // try a different seed value if default guess fails
    if (!solution_found)
    {
      // try twice the grid maxima
      //result = this_eos->get_tbqs_maxima_no_ext();
      //std::for_each(result.begin(), result.end(), [](double &c){ c *= 2.0; });
      const double hc = constants::hbarc_MeVfm;
      result = std::vector<double>({5001.0/hc,0.0,0.0,0.0});
      std::cout << "     - seed: "
                << result[0]*hc << "   " << result[1]*hc << "   "
                << result[2]*hc << "   " << result[3]*hc << std::endl;
      solution_found
        = rootfinder.find_eBSQ_root( ein, Bin, Sin, Qin, this_eos->eBSQ,
                                     this_eos->tbqs_minima, this_eos->tbqs_maxima,
                                     result );
    }

    // stop iterating through available EoSs when solution found
    if (solution_found)
    {
      // any time we update the EoS pointer, we need to specify WHICH EoS we are updating!
      std::cout << " --> found a solution with " << this_eos->name << " EoS!" << std::endl;
      current_eos_name = this_eos->name;
      tbqs( result, this_eos ); // set thermodynamics using solution
      break;
    }
  }

  if (!solution_found)
  {
    std::cout << "No solution found!" << std::endl;
    std::cerr << "No solution found!" << std::endl;
    exit(101);
  }

  // this is set in most recent call to tbqs()
  return entrVal;

}
////////////////////////////////////////////////////////////////////////////////

