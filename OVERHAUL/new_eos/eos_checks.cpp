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
#include "toy_thermo.h"

using namespace constants;

using std::vector;
using std::string;

void EquationOfState::check_EoS_derivatives()
{
	// Check various complicated EoS derivatives assuming simple equation of state
	// Compare with numerical checks from Mathematica to ensure error-free execution

	double Tcheck = toy_thermo::toy_T_scale/2.0;
	double muBcheck = toy_thermo::toy_muB_scale/3.0;
	double muScheck = toy_thermo::toy_muS_scale/7.0;
	double muQcheck = toy_thermo::toy_muQ_scale/5.0;
	/*for (double Tcheck = 30.0; Tcheck <= 800.001; Tcheck += 5.0)
	for (double muBcheck = -450.0; muBcheck <= 450.001; muBcheck += 10.0)
	for (double muScheck = -450.0; muScheck <= 450.001; muScheck += 10.0)
	for (double muQcheck = -450.0; muQcheck <= 450.001; muQcheck += 10.0)*/
	{
		// reset (T,muB,muQ,muS) coordinates
		tbqs(Tcheck/hbarc_MeVfm, muBcheck/hbarc_MeVfm, muQcheck/hbarc_MeVfm, muScheck/hbarc_MeVfm);	// note order of Q and S!
		cout << Tcheck << "   " << muBcheck << "   " << muScheck << "   " << muQcheck << "\n\t"
				<< p() << "   " << s() << "   " << e() << "   "
				<< B() << "   " << S() << "   " << Q() << "\n\t\t"
				<< dt2 << "   " << dtdb << "   " << dtds << "   " << dtdq << "\n\t\t"
				<< dtdb << "   " << db2 << "   " << dbds << "   " << dbdq << "\n\t\t"
				<< dtds << "   " << dbds << "   " << ds2 << "   " << dsdq << "\n\t\t"
				<< dtdq << "   " << dbdq << "   " << dsdq << "   " << dq2 << "\n\t"
				//<< dentr_dt() << "   " << dentr_dmub() << "   "
				//<< dentr_dmus() << "   " << dentr_dmuq() << "\n\t"
				<< dwds() << "   " << dwdB() << "   " << dwdS() << "   " << dwdQ() << "\n";
	}

	return;
}



void EquationOfState::get_toy_thermo(double point[], double thermodynamics[])
{
	const double Tsol = point[0], muBsol = point[1], muSsol = point[2], muQsol = point[3];
	double POut = toy_thermo::p(Tsol, muBsol, muQsol, muSsol);
	double sOut = toy_thermo::s(Tsol, muBsol, muQsol, muSsol);
	double BOut = toy_thermo::B(Tsol, muBsol, muQsol, muSsol);
	double SOut = toy_thermo::S(Tsol, muBsol, muQsol, muSsol);
	double QOut = toy_thermo::Q(Tsol, muBsol, muQsol, muSsol);
	double eOut = (sOut*Tsol + muBsol*BOut + muQsol*QOut + muSsol*SOut)/hbarc_MeVfm - POut;


	//Thermodynamics
	thermodynamics[0]  = POut;
	thermodynamics[1]  = sOut;
	thermodynamics[2]  = BOut;
	thermodynamics[3]  = SOut;
	thermodynamics[4]  = QOut;
	thermodynamics[5]  = eOut;
	thermodynamics[6]  = 0.0;	// not currently checking this
				
	//Second Order Derivatives (prefactor converts to physical susceptibilities in fm^-2)
	thermodynamics[7]  = toy_thermo::P2B2(Tsol, muBsol, muQsol, muSsol);
	thermodynamics[8]  = toy_thermo::P2Q2(Tsol, muBsol, muQsol, muSsol);
	thermodynamics[9]  = toy_thermo::P2S2(Tsol, muBsol, muQsol, muSsol);
	
	thermodynamics[10] = toy_thermo::P2BQ(Tsol, muBsol, muQsol, muSsol);
	thermodynamics[11] = toy_thermo::P2BS(Tsol, muBsol, muQsol, muSsol);
	thermodynamics[12] = toy_thermo::P2QS(Tsol, muBsol, muQsol, muSsol);
	
	thermodynamics[13] = toy_thermo::P2TB(Tsol, muBsol, muQsol, muSsol);
	thermodynamics[14] = toy_thermo::P2TQ(Tsol, muBsol, muQsol, muSsol);
	thermodynamics[15] = toy_thermo::P2TS(Tsol, muBsol, muQsol, muSsol);
	thermodynamics[16] = toy_thermo::P2T2(Tsol, muBsol, muQsol, muSsol);
}

