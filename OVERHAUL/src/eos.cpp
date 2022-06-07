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
#include "../include/eos_delaunay.h"
#include "../src/eos_derivatives.cpp"
#include "../src/eos_initialization.cpp"

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

//cout << "before" << endl;
//cout << "THERMO DUMP: " << pVal << "   " << entrVal << "   " << BVal << "   "
//      << SVal << "   " << QVal << "   " << eVal << "   " << cs2Val << "   "
//      << db2 << "   " << dq2 << "   " << ds2 << "   " << dbdq << "   "
//      << dbds << "   " << dsdq << "   " << dtdb << "   " << dtdq << "   "
//      << dtds << "   " << dt2 << "   " << peos->name << endl;
//cout << "after" << endl;

/*if (true) exit(1);*/
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
//  std::cout << __PRETTY_FUNCTION__ << "::" << __LINE__ << ": "
//            << "using " << eos_name << "!" << std::endl;
  tbqs(T_freeze_out_at_mu_eq_0, 0, 0, 0, eos_name);
//cout << __PRETTY_FUNCTION__ << "::" << __LINE__ << endl;
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
  bool solution_found = false;
  tbqs( result, "table" );
  return solution_found;
}


////////////////////////////////////////////////
bool EquationOfState::find_root_with_seed(
      const string & e_or_s_mode, double e_or_s_in,
      double Bin, double Sin, double Qin,
      pEoS_base this_eos, vector<double> & result )
{
    const double hc = constants::hbarc_MeVfm;
    if ( VERBOSE > 2 )
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
bool EquationOfState::rootfinder_update_s(double sin, double Bin, double Sin, double Qin)
{
  const double hc = constants::hbarc_MeVfm;

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
    solution_found = find_root_with_seed( "entropy", sin, Bin, Sin, Qin, this_eos, result );


    ////////////////////////////////////////////////
    // if rootfinder fails, try a different seed
    if (!solution_found)
    {
      // try default seed at zero density
      result = vector<double>({tbqsPosition[0],0.0,0.0,0.0});
      solution_found = find_root_with_seed( "entropy", sin, Bin, Sin, Qin, this_eos, result );
    }


    ////////////////////////////////////////////////
    // if rootfinder fails, try a different seed
    if (!solution_found)
    {
      // try forced seed
      result = vector<double>({1000.0/hc,0.0,0.0,0.0});
      solution_found = find_root_with_seed( "entropy", sin, Bin, Sin, Qin, this_eos, result );
    }


    ////////////////////////////////////////////////
    // if rootfinder fails, try another different seed
    if (!solution_found)
    {
      // try conformal diagonal seed
      auto conformal_diagonal_EoS = std::dynamic_pointer_cast<EoS_conformal_diagonal>
                                        ( chosen_EOS_map["conformal_diagonal"] );

      result = conformal_diagonal_EoS->get_tbqs_seed_from_sBSQ
                                        ( sin, Bin, Sin, Qin );

      solution_found = find_root_with_seed( "entropy", sin, Bin, Sin, Qin,
                                            this_eos, result );
    }


    ////////////////////////////////////////////////
    // if all else fails and relevant flag is set,
    // ZERO CHARGE DENSITIES
    if (zero_unsolvable_charge_densities && !solution_found)
    {
      result = vector<double>({1000.0/hc,0.0,0.0,0.0});
      solution_found = find_root_with_seed( "entropy", sin, 0.0, 0.0, 0.0,
                                            this_eos, result );
    }



    ////////////////////////////////////////////////
    // stop iterating through available EoSs when solution found
    if (solution_found)
    {
      tbqs( result, this_eos ); // set thermodynamics using solution

      if ( VERBOSE > 2 )
      {
        std::cout << " --> found a solution with "
                  << this_eos->name << " EoS!" << std::endl;
        std::cout << " --> solution has (s,B,S,Q) = "
                  << s() << "   " << B() << "   " << S() << "   " << Q()
                  << std::endl;
      }
      break;
    }
  }

  if (!solution_found)
  {
    std::cout << "No solution found!" << std::endl;
    std::cerr << "No solution found!" << std::endl;
    std::cout << "Last attempted EoS: " << eos_currently_trying << endl;
    std::cout << "Failed to find a solution for (s,B,S,Q) = "
              << sin << "   " << Bin << "   " << Sin << "   " << Qin << std::endl;
    exit(101);
  }

  return solution_found;
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
  const double hc = constants::hbarc_MeVfm;

  // used for seed value in rootfinder
  vector<double> result;
  solution_found = false;

  ///////////////////////////////////////////////////////////
  // try each EoS in turn
  for ( const auto & this_eos : chosen_EOSs )
  {

    /////////////////////////////////////////////////////////
    // try forced seed first
    result = vector<double>({1000.0/hc,0.0,0.0,0.0});
    solution_found = find_root_with_seed( "energy", ein, Bin, Sin, Qin,
                                          this_eos, result );


    /////////////////////////////////////////////////////////
    // try conformal diagonal seed next
    if (!solution_found)
    {
      // try conformal diagonal seed
      auto conformal_diagonal_EoS = std::dynamic_pointer_cast<EoS_conformal_diagonal>
                                        ( chosen_EOS_map["conformal_diagonal"] );

      result = conformal_diagonal_EoS->get_tbqs_seed_from_eBSQ
                                        ( ein, Bin, Sin, Qin );

      solution_found = find_root_with_seed( "energy", ein, Bin, Sin, Qin,
                                            this_eos, result );
    }

    ////////////////////////////////////////////////
    // if all else fails and relevant flag is set,
    // ZERO CHARGE DENSITIES
    if (zero_unsolvable_charge_densities && !solution_found)
    {
      result = vector<double>({1000.0/hc,0.0,0.0,0.0});
      solution_found = find_root_with_seed( "energy", ein, 0.0, 0.0, 0.0,
                                            this_eos, result );
    }


    /////////////////////////////////////////////////////////
    // stop iterating through available EoSs when solution found
    if (solution_found)
    {
      tbqs( result, this_eos ); // set thermodynamics using solution

      if ( VERBOSE > 2 )
      {
        std::cout << " --> found a solution with "
                  << this_eos->name << " EoS!" << std::endl;
        std::cout << " --> solution has (e,B,S,Q) = "
                  << e() << "   " << B() << "   " << S() << "   " << Q() << endl;
      }
      break;
    }
  }

  if (!solution_found)
  {
    std::cout << "No solution found!" << std::endl;
    std::cerr << "No solution found!" << std::endl;
    std::cout << "Failed to find a solution for (e,B,S,Q) = "
              << ein << "   " << Bin << "   " << Sin << "   " << Qin << std::endl;
    exit(101);
  }

  // this is set in most recent call to tbqs()
  return entrVal;

}
////////////////////////////////////////////////////////////////////////////////

