#include "eos.h"

//#include "read_in_hdf/read_in_hdf.h"
#include "Stopwatch.h"
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

//EoS constructor
EquationOfState::EquationOfState(string quantityFile, string derivFile)
{
    init(quantityFile, derivFile);
}

//EoS default constructor. This function exists to satisfy the compiler
//This function should never be called unless init is called directly afterward
EquationOfState::EquationOfState() {}

void EquationOfState::init()
{
  cout << "Attempting read in of EoS from "
        << quantity_file << " and " << deriv_file << endl;
  init( quantity_file, deriv_file );
}

void EquationOfState::init(string quantityFile, string derivFile)
{
	tbqsPosition.resize(4);

	/*if ( check_derivatives )
	{
		cout << "Running EoS in test mode; checking derivatives!" << endl;
		check_EoS_derivatives();
		cout << "All tests completed!  Exiting." << endl;
		exit(-1);
	}*/

	cout << "Initializing EoS C library" << endl;
	initialize("/projects/jnorhos/BSQ/EoS_BQS_Derivatives/Coefficients_Parameters.dat");

	std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
	init_grid_ranges_only(quantityFile, derivFile);

	cout << "Initialize Delaunay interpolators" << endl;
	e_delaunay.init(quantityFile, 0);		// 0 - energy density
	entr_delaunay.init(quantityFile, 1);	// 1 - entropy density

	return;
}

void EquationOfState::init_grid_ranges_only(string quantityFile, string derivFile)
{
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    std::ifstream dataFile;
    dataFile.open(quantityFile);

    double tit, muBit, muQit, muSit, pit, entrit, bit, sit, qit, eit, cs2it;

    int count = 0;
    double hc = hbarc_MeVfm;
    while (dataFile >> tit >> muBit >> muQit >> muSit
			>> pit >> entrit >> bit >> sit >> qit
			>> eit >> cs2it)
    {

		// Christopher Plumberg:
		// put T and mu_i in units of 1/fm
		tit   /= hc;
		muBit /= hc;
		muSit /= hc;
		muQit /= hc;

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

    dataFile.close();

	std::cout << "All initializations finished!" << std::endl;

    return;
}

void EquationOfState::tbqs(double setT, double setmuB, double setmuQ, double setmuS)
{
	//if ( !check_derivatives )
	//{
		if(setT < minT || setT > maxT) {
			std::cout << "T = " << setT << " is out of range. Valid values are between ["
				<< minT << "," << maxT << "]" << std::endl;
			return;
		}
		if(setmuB < minMuB || setmuB > maxMuB) {
			std::cout << "muB = " << setmuB << " is out of range. Valid values are between ["
				<< minMuB << "," << maxMuB << "]" << std::endl;
			return;
		}
		if(setmuQ < minMuQ || setmuQ > maxMuQ) {
			std::cout << "muQ = " << setmuQ << " is out of range. Valid values are between ["
				<< minMuQ << "," << maxMuQ << "]" << std::endl;
			return;
		}
		if(setmuS < minMuS || setmuS > maxMuS) {
			std::cout << "muS = " << setmuS << " is out of range. Valid values are between ["
				<< minMuS << "," << maxMuS << "]" << std::endl;
			return;
		}
	//}
	tbqsPosition[0] = setT;
	tbqsPosition[1] = setmuB;
	tbqsPosition[2] = setmuQ;
	tbqsPosition[3] = setmuS;

	// EXPECTS UNITS OF MEV!!!
	double phase_diagram_point[4]	// NOTE: S <<-->> Q swapped!!!
			= {setT*hbarc_MeVfm, setmuB*hbarc_MeVfm, setmuS*hbarc_MeVfm, setmuQ*hbarc_MeVfm};
	double thermodynamics[17];
	//if ( check_derivatives )
	//	get_toy_thermo(phase_diagram_point, thermodynamics);
	//else
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




double EquationOfState::delaunay_s_out(double ein, double Bin, double Sin, double Qin)
{
  
}




double EquationOfState::rootfinder_s_out(double ein, double Bin, double Sin, double Qin)
{
  
}
