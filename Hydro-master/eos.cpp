#include "eos.h"

#include "read_in_hdf/read_in_hdf.h"
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


using std::vector;
using std::string;

// Compile:         gcc eos4D.cpp -c -I /usr/include/eigen3 -Lsplinter/build -lm -lgsl -lgslcblas -lstdc++ -lsplinter-3-0

constexpr bool use_exact = true;
constexpr bool accept_nearest_neighbor = false;
constexpr bool discard_unsolvable_charge_densities = false;
constexpr bool check_derivatives = true;

constexpr size_t STEPS = 100000;
constexpr int VERBOSE = 0;
constexpr double TOLERANCE = 1e-12;

namespace toy_thermo
{
	const double hc = 197.3;
	const double toy_T_scale = 150.0, toy_muB_scale = 101.0, toy_muS_scale = 103.0, toy_muQ_scale = 113.0;

	double p(double T, double muB, double muQ, double muS)
	{
		double x = (T/toy_T_scale)*(T/toy_T_scale) + (muB/toy_muB_scale)*(muB/toy_muB_scale)
					 + (muS/toy_muS_scale)*(muS/toy_muS_scale) + (muQ/toy_muQ_scale)*(muQ/toy_muQ_scale);
		return x*x;
	}
	
	double s(double T, double muB, double muQ, double muS)
	{
		double x = (T/toy_T_scale)*(T/toy_T_scale) + (muB/toy_muB_scale)*(muB/toy_muB_scale)
					 + (muS/toy_muS_scale)*(muS/toy_muS_scale) + (muQ/toy_muQ_scale)*(muQ/toy_muQ_scale);
		return 4.0*hc*T*x/(toy_T_scale*toy_T_scale);
	}

	double B(double T, double muB, double muQ, double muS)
	{
		double x = (T/toy_T_scale)*(T/toy_T_scale) + (muB/toy_muB_scale)*(muB/toy_muB_scale)
					 + (muS/toy_muS_scale)*(muS/toy_muS_scale) + (muQ/toy_muQ_scale)*(muQ/toy_muQ_scale);
		return 4.0*hc*muB*x/(toy_muB_scale*toy_muB_scale);
	}

	double S(double T, double muB, double muQ, double muS)
	{
		double x = (T/toy_T_scale)*(T/toy_T_scale) + (muB/toy_muB_scale)*(muB/toy_muB_scale)
					 + (muS/toy_muS_scale)*(muS/toy_muS_scale) + (muQ/toy_muQ_scale)*(muQ/toy_muQ_scale);
		return 4.0*hc*muS*x/(toy_muS_scale*toy_muS_scale);
	}

	double Q(double T, double muB, double muQ, double muS)
	{
		double x = (T/toy_T_scale)*(T/toy_T_scale) + (muB/toy_muB_scale)*(muB/toy_muB_scale)
					 + (muS/toy_muS_scale)*(muS/toy_muS_scale) + (muQ/toy_muQ_scale)*(muQ/toy_muQ_scale);
		return 4.0*hc*muQ*x/(toy_muQ_scale*toy_muQ_scale);
	}

	double P2B2(double T, double muB, double muQ, double muS)
	{
		double x = (T/toy_T_scale)*(T/toy_T_scale) + (muB/toy_muB_scale)*(muB/toy_muB_scale)
					 + (muS/toy_muS_scale)*(muS/toy_muS_scale) + (muQ/toy_muQ_scale)*(muQ/toy_muQ_scale);
		return (8.0*muB*muB/(toy_muB_scale*toy_muB_scale*toy_muB_scale*toy_muB_scale)
				+ 4.0*x/(toy_muB_scale*toy_muB_scale))*hc*hc;
	}
	double P2Q2(double T, double muB, double muQ, double muS)
	{
		double x = (T/toy_T_scale)*(T/toy_T_scale) + (muB/toy_muB_scale)*(muB/toy_muB_scale)
					 + (muS/toy_muS_scale)*(muS/toy_muS_scale) + (muQ/toy_muQ_scale)*(muQ/toy_muQ_scale);
		return (8.0*muQ*muQ/(toy_muQ_scale*toy_muQ_scale*toy_muQ_scale*toy_muQ_scale)
				+ 4.0*x/(toy_muQ_scale*toy_muQ_scale))*hc*hc;
	}
	double P2S2(double T, double muB, double muQ, double muS)
	{
		double x = (T/toy_T_scale)*(T/toy_T_scale) + (muB/toy_muB_scale)*(muB/toy_muB_scale)
					 + (muS/toy_muS_scale)*(muS/toy_muS_scale) + (muQ/toy_muQ_scale)*(muQ/toy_muQ_scale);
		return (8.0*muS*muS/(toy_muS_scale*toy_muS_scale*toy_muS_scale*toy_muS_scale)
				+ 4.0*x/(toy_muS_scale*toy_muS_scale))*hc*hc;
	}
	
	double P2BQ(double T, double muB, double muQ, double muS)
			{ return (8.0*muB*muQ/(toy_muB_scale*toy_muB_scale*toy_muQ_scale*toy_muQ_scale))*hc*hc; }
	double P2BS(double T, double muB, double muQ, double muS)
			{ return (8.0*muB*muS/(toy_muB_scale*toy_muB_scale*toy_muS_scale*toy_muS_scale))*hc*hc; }
	double P2QS(double T, double muB, double muQ, double muS)
			{ return (8.0*muS*muQ/(toy_muS_scale*toy_muS_scale*toy_muQ_scale*toy_muQ_scale))*hc*hc; }
	
	double P2TB(double T, double muB, double muQ, double muS)
			{ return (8.0*T*muB/(toy_T_scale*toy_T_scale*toy_muB_scale*toy_muB_scale))*hc*hc; }
	double P2TQ(double T, double muB, double muQ, double muS)
			{ return (8.0*T*muQ/(toy_T_scale*toy_T_scale*toy_muQ_scale*toy_muQ_scale))*hc*hc; }
	double P2TS(double T, double muB, double muQ, double muS)
			{ return (8.0*T*muS/(toy_T_scale*toy_T_scale*toy_muS_scale*toy_muS_scale))*hc*hc; }
	double P2T2(double T, double muB, double muQ, double muS)
	{
		double x = (T/toy_T_scale)*(T/toy_T_scale) + (muB/toy_muB_scale)*(muB/toy_muB_scale)
					 + (muS/toy_muS_scale)*(muS/toy_muS_scale) + (muQ/toy_muQ_scale)*(muQ/toy_muQ_scale);
		return (8.0*T*T/(toy_T_scale*toy_T_scale*toy_T_scale*toy_T_scale)
				+ 4.0*x/(toy_T_scale*toy_T_scale))*hc*hc;
	}



}

//EoS constructor
eos::eos(string quantityFile, string derivFile)
{
    init(quantityFile, derivFile);
}

//EoS default constructor. This function exists to satisfy the compiler
//This function should never be called unless init is called directly afterward
eos::eos() {}

void eos::init(string quantityFile, string derivFile)
{
	tbqsPosition.resize(4);

	if ( check_derivatives )
	{
		cout << "Running EoS in test mode; checking derivatives!" << endl;
		check_EoS_derivatives();
		cout << "All tests completed!  Exiting." << endl;
		exit(-1);
	}

	cout << "Initializing EoS C library" << endl;
	initialize("/projects/jnorhos/BSQ/EoS_BQS_Derivatives/Coefficients_Parameters.dat");

	std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
	init_grid_ranges_only(quantityFile, derivFile);

	cout << "Initialize Delaunay interpolators" << endl;
	e_delaunay.init(quantityFile, 0);		// 0 - energy density
	entr_delaunay.init(quantityFile, 1);	// 1 - entropy density

	return;
}

void eos::check_EoS_derivatives()
{
	// Check various complicated EoS derivatives assuming simple equation of state
	// Compare with numerical checks from Mathematica to ensure error-free execution

	double Tcheck = 0.5*toy_thermo::toy_T_scale;
	double muBcheck = 0.0*toy_thermo::toy_muB_scale;
	double muScheck = 0.0*toy_thermo::toy_muS_scale;
	double muQcheck = 0.0*toy_thermo::toy_muQ_scale;
	/*for (double Tcheck = 30.0; Tcheck <= 800.001; Tcheck += 5.0)
	for (double muBcheck = -450.0; muBcheck <= 450.001; muBcheck += 10.0)
	for (double muScheck = -450.0; muScheck <= 450.001; muScheck += 10.0)
	for (double muQcheck = -450.0; muQcheck <= 450.001; muQcheck += 10.0)*/
	{
		// reset (T,muB,muQ,muS) coordinates
		tbqs(Tcheck/197.3, muBcheck/197.3, muQcheck/197.3, muScheck/197.3);	// note order of Q and S!
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



void eos::get_toy_thermo(double point[], double thermodynamics[])
{
	const double Tsol = point[0], muBsol = point[1], muSsol = point[2], muQsol = point[3];
	double POut = toy_thermo::p(Tsol, muBsol, muQsol, muSsol);
	double sOut = toy_thermo::s(Tsol, muBsol, muQsol, muSsol);
	double BOut = toy_thermo::B(Tsol, muBsol, muQsol, muSsol);
	double SOut = toy_thermo::S(Tsol, muBsol, muQsol, muSsol);
	double QOut = toy_thermo::Q(Tsol, muBsol, muQsol, muSsol);
	double eOut = (sOut*Tsol + muBsol*BOut + muQsol*QOut + muSsol*SOut)/197.3 - POut;


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



void eos::init_grid_ranges_only(string quantityFile, string derivFile)
{
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    std::ifstream dataFile;
    dataFile.open(quantityFile);

    double tit, muBit, muQit, muSit, pit, entrit, bit, sit, qit, eit, cs2it;

    int count = 0;
    double hbarc = 197.3;
    while (dataFile >> tit >> muBit >> muQit >> muSit
			>> pit >> entrit >> bit >> sit >> qit
			>> eit >> cs2it)
    {

		// Christopher Plumberg:
		// put T and mu_i in units of 1/fm
		tit   /= hbarc;
		muBit /= hbarc;
		muSit /= hbarc;
		muQit /= hbarc;

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

void eos::tbqs(double setT, double setmuB, double setmuQ, double setmuS)
{
	if ( !check_derivatives )
	{
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
	}
	tbqsPosition[0] = setT;
	tbqsPosition[1] = setmuB;
	tbqsPosition[2] = setmuQ;
	tbqsPosition[3] = setmuS;

	// EXPECTS UNITS OF MEV!!!
	double phase_diagram_point[4]	// NOTE: S <<-->> Q swapped!!!
			= {setT*197.3, setmuB*197.3, setmuS*197.3, setmuQ*197.3};
	double thermodynamics[17];
	if ( check_derivatives )
		get_toy_thermo(phase_diagram_point, thermodynamics);
	else
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


double eos::T() { return tbqsPosition[0]; }
double eos::muB() { return tbqsPosition[1]; }
double eos::muQ() { return tbqsPosition[2]; }
double eos::muS() { return tbqsPosition[3]; }

double eos::p() { return pVal; }
double eos::s() { return entrVal; }
double eos::B() { return BVal; }
double eos::S() { return SVal; }
double eos::Q() { return QVal; }
double eos::e() { return eVal; }
double eos::cs2() { return cs2Val; }
double eos::w() { return eVal + pVal; }


double eos::dwds()
{
	double charge_terms	/*if charge densities are not all zero*/
			= ( abs(BVal)>1e-10 || abs(SVal)>1e-10 || abs(QVal)>1e-10 ) ?
			  BVal/dentr_dmub() + QVal/dentr_dmuq() + SVal/dentr_dmus() : 0.0;

	cout << endl << endl << "inside dwds(): "
		<< T() << "   " << entrVal << "   " << dentr_dt() << "   " << charge_terms << endl << endl;

    return T() + entrVal/dentr_dt() + charge_terms;
}

double eos::dwdB()
{
	double charge_terms	/*if charge densities are not all zero*/
			= ( abs(BVal)>1e-10 || abs(SVal)>1e-10 || abs(QVal)>1e-10 ) ?
			  entrVal/db_dt() + BVal/db_dmub() + QVal/db_dmuq() + SVal/db_dmus() : 0.0;

    return muB() + charge_terms;
}

double eos::dwdS()
{
	double charge_terms	/*if charge densities are not all zero*/
			= ( abs(BVal)>1e-10 || abs(SVal)>1e-10 || abs(QVal)>1e-10 ) ?
			  entrVal/ds_dt() + BVal/ds_dmub() + QVal/ds_dmuq() + SVal/ds_dmus() : 0.0;

    return muS() + charge_terms;
}

double eos::dwdQ()
{
	double charge_terms	/*if charge densities are not all zero*/
			= ( abs(BVal)>1e-10 || abs(SVal)>1e-10 || abs(QVal)>1e-10 ) ?
			  entrVal/dq_dt() + BVal/dq_dmub() + QVal/dq_dmuq() + SVal/dq_dmus() : 0.0;

    return muQ() + charge_terms;
}

double eos::dentr_dt()   { return calc_term_1();        }
double eos::dentr_dmub() { return calc_term_2("b");     }
double eos::dentr_dmuq() { return calc_term_2("q");     }
double eos::dentr_dmus() { return calc_term_2("s");     }
double eos::db_dt()      { return calc_term_3("b");     }
double eos::db_dmub()    { return calc_term_4("b","b"); }
double eos::db_dmuq()    { return calc_term_4("b","q"); }
double eos::db_dmus()    { return calc_term_4("b","s"); }
double eos::ds_dt()      { return calc_term_3("s");     }
double eos::ds_dmub()    { return calc_term_4("s","b"); }
double eos::ds_dmuq()    { return calc_term_4("s","q"); }
double eos::ds_dmus()    { return calc_term_4("s","s"); }
double eos::dq_dt()      { return calc_term_3("q");     }
double eos::dq_dmub()    { return calc_term_4("q","b"); }
double eos::dq_dmuq()    { return calc_term_4("q","q"); }
double eos::dq_dmus()    { return calc_term_4("q","s"); }

double eos::calc_term_1() {
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    gsl_vector *v = gsl_vector_alloc(3);
    gsl_matrix *m = gsl_matrix_alloc(3,3);

    gsl_vector_set(v,0,dtdb);
    gsl_vector_set(v,1,dtds);
    gsl_vector_set(v,2,dtdq);

/*if (verbose)
{
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
	cout << "calc_term_1 check: v = ";
	for (int iGSL = 0; iGSL < 3; iGSL++)
		cout << gsl_vector_set(v,iGSL);
	cout << endl;
}*/

    gsl_matrix_set(m,0,0,db2);
    gsl_matrix_set(m,0,1,dbds);
    gsl_matrix_set(m,0,2,dbdq);
    gsl_matrix_set(m,1,0,dbds);
    gsl_matrix_set(m,1,1,ds2);
    gsl_matrix_set(m,1,2,dsdq);
    gsl_matrix_set(m,2,0,dbdq);
    gsl_matrix_set(m,2,1,dsdq);
    gsl_matrix_set(m,2,2,dq2);

    double toReturn = dt2 - deriv_mult_aTm_1b(v,m,v);

/*if (verbose)
{
	cout << "calc_term_1 check: m = ";
	for (int iGSL = 0; iGSL < 3; iGSL++)
	{
		for (int jGSL = 0; jGSL < 3; jGSL++)
			cout << "   " << gsl_matrix_set(m,iGSL,jGSL);
		cout << endl;
	}
	cout << "calc_term_1 check: " << toReturn << "   "
			<< dt2 << "   " << deriv_mult_aTm_1b(v,m,v) << endl;
	cout << "----------------------------------------"
			"----------------------------------------" << endl;
}*/

    gsl_matrix_free(m);
    gsl_vector_free(v);
    return toReturn;
}

double eos::calc_term_2(string i_char) {
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__
								 << ": i_char = " << i_char << std::endl;
    gsl_vector *a = gsl_vector_alloc(3);
    gsl_matrix *m = gsl_matrix_alloc(3,3);
    gsl_vector *b = gsl_vector_alloc(3);
    double toReturn = 0;

    if (i_char == "b") {
        gsl_vector_set(a,0,dt2);
        gsl_vector_set(a,1,dtds);
        gsl_vector_set(a,2,dtdq);

        gsl_vector_set(b,0,db2);
        gsl_vector_set(b,1,dbds);
        gsl_vector_set(b,2,dbdq);

        gsl_matrix_set(m,0,0,dtdb);
        gsl_matrix_set(m,0,1,dbds);
        gsl_matrix_set(m,0,2,dbdq);
        gsl_matrix_set(m,1,0,dtds);
        gsl_matrix_set(m,1,1,ds2);
        gsl_matrix_set(m,1,2,dsdq);
        gsl_matrix_set(m,2,0,dtdq);
        gsl_matrix_set(m,2,1,dsdq);
        gsl_matrix_set(m,2,2,dq2);

        toReturn = dtdb - deriv_mult_aTm_1b(a,m,b);
    } else if (i_char == "s") {
        gsl_vector_set(a,0,dt2);
        gsl_vector_set(a,1,dtdb);
        gsl_vector_set(a,2,dtdq);

        gsl_vector_set(b,0,dbds);
        gsl_vector_set(b,1,ds2);
        gsl_vector_set(b,2,dsdq);

        gsl_matrix_set(m,0,0,dtdb);
        gsl_matrix_set(m,0,1,db2);
        gsl_matrix_set(m,0,2,dbdq);
        gsl_matrix_set(m,1,0,dtds);
        gsl_matrix_set(m,1,1,dbds);
        gsl_matrix_set(m,1,2,dsdq);
        gsl_matrix_set(m,2,0,dtdq);
        gsl_matrix_set(m,2,1,dbdq);
        gsl_matrix_set(m,2,2,dq2);

        toReturn = dtds - deriv_mult_aTm_1b(a,m,b);
    } else if (i_char == "q") {
        gsl_vector_set(a,0,dt2);
        gsl_vector_set(a,1,dtdb);
        gsl_vector_set(a,2,dtds);

        gsl_vector_set(b,0,dbdq);
        gsl_vector_set(b,1,dsdq);
        gsl_vector_set(b,2,dq2);

        gsl_matrix_set(m,0,0,dtdb);
        gsl_matrix_set(m,0,1,db2);
        gsl_matrix_set(m,0,2,dbdq);
        gsl_matrix_set(m,1,0,dtds);
        gsl_matrix_set(m,1,1,dbds);
        gsl_matrix_set(m,1,2,dsdq);
        gsl_matrix_set(m,2,0,dtdq);
        gsl_matrix_set(m,2,1,dbdq);
        gsl_matrix_set(m,2,2,dq2);

        toReturn = dtdq - deriv_mult_aTm_1b(a,m,b);
    } else {
        std::cout << "Error calculating derivative term 2" << std::endl;
    }

    gsl_vector_free(a);
    gsl_matrix_free(m);
    gsl_vector_free(b);
    return toReturn;
}

double eos::calc_term_3(string i_char) {
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__
								 << ": i_char = " << i_char << std::endl;
    gsl_vector *a = gsl_vector_alloc(3);
    gsl_matrix *m = gsl_matrix_alloc(3,3);
    gsl_vector *b = gsl_vector_alloc(3);
    double toReturn = 0;

    if (i_char == "b") {
        gsl_vector_set(a,0,db2);
        gsl_vector_set(a,1,dbds);
        gsl_vector_set(a,2,dbdq);

        gsl_vector_set(b,0,dt2);
        gsl_vector_set(b,1,dtds);
        gsl_vector_set(b,2,dtdq);

        gsl_matrix_set(m,0,0,dtdb);
        gsl_matrix_set(m,0,1,dtds);
        gsl_matrix_set(m,0,2,dtdq);
        gsl_matrix_set(m,1,0,dbds);
        gsl_matrix_set(m,1,1,ds2);
        gsl_matrix_set(m,1,2,dsdq);
        gsl_matrix_set(m,2,0,dbdq);
        gsl_matrix_set(m,2,1,dsdq);
        gsl_matrix_set(m,2,2,dq2);

        toReturn = dtdb - deriv_mult_aTm_1b(a,m,b);
    } else if (i_char == "s") {
        gsl_vector_set(a,0,dbds);
        gsl_vector_set(a,1,ds2);
        gsl_vector_set(a,2,dsdq);

        gsl_vector_set(b,0,dt2);
        gsl_vector_set(b,1,dtdb);
        gsl_vector_set(b,2,dtdq);

        gsl_matrix_set(m,0,0,dtdb);
        gsl_matrix_set(m,0,1,dtds);
        gsl_matrix_set(m,0,2,dtdq);
        gsl_matrix_set(m,1,0,db2);
        gsl_matrix_set(m,1,1,dbds);
        gsl_matrix_set(m,1,2,dbdq);
        gsl_matrix_set(m,2,0,dbdq);
        gsl_matrix_set(m,2,1,dsdq);
        gsl_matrix_set(m,2,2,dq2);

        toReturn = dtds - deriv_mult_aTm_1b(a,m,b);
    } else if (i_char == "q") {
        gsl_vector_set(a,0,dbdq);
        gsl_vector_set(a,1,dsdq);
        gsl_vector_set(a,2,dq2);

        gsl_vector_set(b,0,dt2);
        gsl_vector_set(b,1,dtdb);
        gsl_vector_set(b,2,dtds);

        gsl_matrix_set(m,0,0,dtdb);
        gsl_matrix_set(m,0,1,dtds);
        gsl_matrix_set(m,0,2,dtdq);
        gsl_matrix_set(m,1,0,db2);
        gsl_matrix_set(m,1,1,dbds);
        gsl_matrix_set(m,1,2,dbdq);
        gsl_matrix_set(m,2,0,dbds);
        gsl_matrix_set(m,2,1,ds2);
        gsl_matrix_set(m,2,2,dsdq);

        toReturn = dtdq - deriv_mult_aTm_1b(a,m,b);
    } else {
        std::cout << "Error calculating derivative term 3" << std::endl;
    }

    gsl_vector_free(a);
    gsl_matrix_free(m);
    gsl_vector_free(b);
    return toReturn;
}

double eos::calc_term_4(string j_char, string i_char) {
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__
								 << ": j_char, i_char = "
								<< j_char << "   " << i_char << std::endl;
    gsl_vector *a = gsl_vector_alloc(3);
    gsl_matrix *m = gsl_matrix_alloc(3,3);
    gsl_vector *b = gsl_vector_alloc(3);
    double toReturn = 0;

    if (i_char == "b") {
        if(j_char == "b") {
            gsl_vector_set(a,0,dtdb);
            gsl_vector_set(a,1,dbds);
            gsl_vector_set(a,2,dbdq);

            gsl_vector_set(b,0,dtdb);
            gsl_vector_set(b,1,dbds);
            gsl_vector_set(b,2,dbdq);

            gsl_matrix_set(m,0,0,dt2);
            gsl_matrix_set(m,0,1,dtds);
            gsl_matrix_set(m,0,2,dtdq);
            gsl_matrix_set(m,1,0,dtds);
            gsl_matrix_set(m,1,1,ds2);
            gsl_matrix_set(m,1,2,dsdq);
            gsl_matrix_set(m,2,0,dtdq);
            gsl_matrix_set(m,2,1,dsdq);
            gsl_matrix_set(m,2,2,dq2);

            toReturn = db2 - deriv_mult_aTm_1b(a,m,b);
        } else if (j_char == "s") {
            gsl_vector_set(a,0,dtds);
            gsl_vector_set(a,1,ds2);
            gsl_vector_set(a,2,dsdq);

            gsl_vector_set(b,0,dtdb);
            gsl_vector_set(b,1,db2);
            gsl_vector_set(b,2,dbdq);

            gsl_matrix_set(m,0,0,dt2);
            gsl_matrix_set(m,0,1,dtds);
            gsl_matrix_set(m,0,2,dtdq);
            gsl_matrix_set(m,1,0,dtdb);
            gsl_matrix_set(m,1,1,dbds);
            gsl_matrix_set(m,1,2,dbdq);
            gsl_matrix_set(m,2,0,dtdq);
            gsl_matrix_set(m,2,1,dsdq);
            gsl_matrix_set(m,2,2,dq2);

            toReturn = dbds - deriv_mult_aTm_1b(a,m,b);
        } else if (j_char == "q") {
            gsl_vector_set(a,0,dtdq);
            gsl_vector_set(a,1,dsdq);
            gsl_vector_set(a,2,dq2);

            gsl_vector_set(b,0,dtdb);
            gsl_vector_set(b,1,db2);
            gsl_vector_set(b,2,dbds);

            gsl_matrix_set(m,0,0,dt2);
            gsl_matrix_set(m,0,1,dtds);
            gsl_matrix_set(m,0,2,dtdq);
            gsl_matrix_set(m,1,0,dtdb);
            gsl_matrix_set(m,1,1,dbds);
            gsl_matrix_set(m,1,2,dbdq);
            gsl_matrix_set(m,2,0,dtds);
            gsl_matrix_set(m,2,1,ds2);
            gsl_matrix_set(m,2,2,dsdq);

            toReturn = dbdq - deriv_mult_aTm_1b(a,m,b);
        } else {
            std::cout << "Error calculating derivative term 4" << std::endl;
        }
    } else if (i_char == "s") {
        if(j_char == "b") {
            gsl_vector_set(a,0,dtdb);
            gsl_vector_set(a,1,db2);
            gsl_vector_set(a,2,dbdq);

            gsl_vector_set(b,0,dtds);
            gsl_vector_set(b,1,ds2);
            gsl_vector_set(b,2,dsdq);

            gsl_matrix_set(m,0,0,dt2);
            gsl_matrix_set(m,0,1,dtdb);
            gsl_matrix_set(m,0,2,dtdq);
            gsl_matrix_set(m,1,0,dtds);
            gsl_matrix_set(m,1,1,dbds);
            gsl_matrix_set(m,1,2,dsdq);
            gsl_matrix_set(m,2,0,dtdq);
            gsl_matrix_set(m,2,1,dbdq);
            gsl_matrix_set(m,2,2,dq2);

            toReturn = dbds - deriv_mult_aTm_1b(a,m,b);
        } else if (j_char == "s") {
            gsl_vector_set(a,0,dtds);
            gsl_vector_set(a,1,dbds);
            gsl_vector_set(a,2,dsdq);

            gsl_vector_set(b,0,dtds);
            gsl_vector_set(b,1,dbds);
            gsl_vector_set(b,2,dsdq);

            gsl_matrix_set(m,0,0,dt2);
            gsl_matrix_set(m,0,1,dtdb);
            gsl_matrix_set(m,0,2,dtdq);
            gsl_matrix_set(m,1,0,dtdb);
            gsl_matrix_set(m,1,1,db2);
            gsl_matrix_set(m,1,2,dbdq);
            gsl_matrix_set(m,2,0,dtdq);
            gsl_matrix_set(m,2,1,dbdq);
            gsl_matrix_set(m,2,2,dq2);

            toReturn = ds2 - deriv_mult_aTm_1b(a,m,b);
        } else if (j_char == "q") {
            gsl_vector_set(a,0,dtdq);
            gsl_vector_set(a,1,dbdq);
            gsl_vector_set(a,2,dq2);

            gsl_vector_set(b,0,dtds);
            gsl_vector_set(b,1,dbds);
            gsl_vector_set(b,2,ds2);

            gsl_matrix_set(m,0,0,dt2);
            gsl_matrix_set(m,0,1,dtdb);
            gsl_matrix_set(m,0,2,dtdq);
            gsl_matrix_set(m,1,0,dtdb);
            gsl_matrix_set(m,1,1,db2);
            gsl_matrix_set(m,1,2,dbdq);
            gsl_matrix_set(m,2,0,dtds);
            gsl_matrix_set(m,2,1,dbds);
            gsl_matrix_set(m,2,2,dsdq);

            toReturn = dsdq - deriv_mult_aTm_1b(a,m,b);
        } else {
            std::cout << "Error calculating derivative term 4" << std::endl;
        }
    } else if (i_char == "q") {
        if(j_char == "b") {
            gsl_vector_set(a,0,dtdb);
            gsl_vector_set(a,1,db2);
            gsl_vector_set(a,2,dbds);

            gsl_vector_set(b,0,dtdq);
            gsl_vector_set(b,1,dsdq);
            gsl_vector_set(b,2,dq2);

            gsl_matrix_set(m,0,0,dt2);
            gsl_matrix_set(m,0,1,dtdb);
            gsl_matrix_set(m,0,2,dtds);
            gsl_matrix_set(m,1,0,dtds);
            gsl_matrix_set(m,1,1,dbds);
            gsl_matrix_set(m,1,2,ds2);
            gsl_matrix_set(m,2,0,dtdq);
            gsl_matrix_set(m,2,1,dbdq);
            gsl_matrix_set(m,2,2,dsdq);

            toReturn = dbdq - deriv_mult_aTm_1b(a,m,b);
        } else if (j_char == "s") {
            gsl_vector_set(a,0,dtds);
            gsl_vector_set(a,1,dbds);
            gsl_vector_set(a,2,ds2);

            gsl_vector_set(b,0,dtdq);
            gsl_vector_set(b,1,dbdq);
            gsl_vector_set(b,2,dq2);

            gsl_matrix_set(m,0,0,dt2);
            gsl_matrix_set(m,0,1,dtdb);
            gsl_matrix_set(m,0,2,dtds);
            gsl_matrix_set(m,1,0,dtdb);
            gsl_matrix_set(m,1,1,db2);
            gsl_matrix_set(m,1,2,dbds);
            gsl_matrix_set(m,2,0,dtdq);
            gsl_matrix_set(m,2,1,dbdq);
            gsl_matrix_set(m,2,2,dsdq);

            toReturn = dsdq - deriv_mult_aTm_1b(a,m,b);
        } else if (j_char == "q") {
            gsl_vector_set(a,0,dtdq);
            gsl_vector_set(a,1,dbdq);
            gsl_vector_set(a,2,dsdq);

            gsl_vector_set(b,0,dtdq);
            gsl_vector_set(b,1,dbdq);
            gsl_vector_set(b,2,dsdq);

            gsl_matrix_set(m,0,0,dt2);
            gsl_matrix_set(m,0,1,dtdb);
            gsl_matrix_set(m,0,2,dtds);
            gsl_matrix_set(m,1,0,dtdb);
            gsl_matrix_set(m,1,1,db2);
            gsl_matrix_set(m,1,2,dbds);
            gsl_matrix_set(m,2,0,dtds);
            gsl_matrix_set(m,2,1,dbds);
            gsl_matrix_set(m,2,2,ds2);

            toReturn = dq2 - deriv_mult_aTm_1b(a,m,b);
        } else {
            std::cout << "Error calculating derivative term 4" << std::endl;
        }
    } else {
        std::cout << "Error calculating derivative term 4" << std::endl;
    }

    gsl_vector_free(a);
    gsl_matrix_free(m);
    gsl_vector_free(b);
    return toReturn;
}

double eos::deriv_mult_aTm_1b(gsl_vector* a, gsl_matrix* m, gsl_vector* b) {
    gsl_permutation *p = gsl_permutation_alloc(3);
    int s;

    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(m, p, &s);

	gsl_set_error_handler_off();

    // Compute the  inverse of the LU decomposition
    gsl_matrix *minv = gsl_matrix_alloc(3, 3);
    int inversion_status = gsl_linalg_LU_invert(m, p, minv);

	if ( inversion_status )	// if an error occurred
	{
		cout << "Current TBQS location: "
				<< 197.3*T() << "   " << 197.3*muB() << "   "
				<< 197.3*muS() << "   " << 197.3*muQ() << endl << endl;

		cout << "Current EoS data:" << endl;
		cout << "pVal = " << pVal << endl
			 << "BVal = " << BVal << endl
			 << "SVal = " << SVal << endl
			 << "QVal = " << QVal << endl
			 << "eVal = " << eVal << endl
			 << "cs2Val = " << cs2Val << endl
			 << "db2 = " << db2 << endl
			 << "ds2 = " << ds2 << endl
			 << "dq2 = " << dq2 << endl
			 << "dt2 = " << dt2 << endl
			 << "dbdq = " << dbdq << endl
			 << "dbds = " << dbds << endl
			 << "dsdq = " << dsdq << endl
			 << "dtdb = " << dtdb << endl
			 << "dtds = " << dtds << endl
			 << "dtdq = " << dtdq << endl
			 << "entrVal = " << entrVal << endl << endl;


		cout << "a=" << endl;
		for (int ii = 0; ii < 3; ii++)
			cout << gsl_vector_get(a, ii) << "   ";
		cout << endl;

		cout << endl << "b=" << endl;
		for (int ii = 0; ii < 3; ii++)
			cout << gsl_vector_get(b, ii) << "   ";
		cout << endl;

		cout << endl << "m=" << endl;
		for (int ii = 0; ii < 3; ii++)
		{
			for (int jj = 0; jj < 3; jj++)
				cout << gsl_matrix_get(m, ii, jj) << "   ";
			cout << endl;
		}
		cout << endl << "minv=" << endl;
		for (int ii = 0; ii < 3; ii++)
		{
			for (int jj = 0; jj < 3; jj++)
				cout << gsl_matrix_get(minv, ii, jj) << "   ";
			cout << endl;
		}
		cout << endl;
		cout << "Exiting!" << endl;
		exit (-1);
	}
	gsl_set_error_handler (NULL);


    gsl_vector *y = gsl_vector_alloc(3);

    // Compute y = m^-1 @ b
    //gsl_blas_dgemv(CblasNoTrans,1,m,b,0,y);
    gsl_blas_dgemv(CblasNoTrans,1,minv,b,0,y);


	/*cout << endl << endl;
	cout << "=============================================" << endl;
	cout << "y=" << endl;
	for (int ii = 0; ii < 3; ii++)
		cout << gsl_vector_get(y, ii) << "   ";
	cout << endl;
	cout << "=============================================" << endl;
	cout << endl << endl;*/


    double toReturn = 0;
    //compute toReturn = aT @ y
    gsl_blas_ddot(a,y,&toReturn);

    gsl_vector_free(y);
    gsl_matrix_free(minv);
    gsl_permutation_free(p);

    return toReturn;
}

double eos::Atable()
{
    Aout=w()-entrVal*dwds();

    return Aout;
}

double eos::cs2out(double Tt) {  //return cs2 given t and mu's=0
    tbqs(Tt, 0.0, 0.0, 0.0);
    return cs2Val;
}

double eos::cs2out(double Tt, double muBin, double muQin, double muSin) {  //return cs2 given t and mu's
    tbqs(Tt, muBin, muQin, muSin);
    return cs2Val;
}

double eos::wfz(double Tt) {   // return e + p for tbqs
    tbqs(Tt, 0.0, 0.0, 0.0);
    return eVal + pVal;
}

double eos::wfz(double Tt, double muBin, double muQin, double muSin) {   // return e + p for tbqs
    tbqs(Tt, muBin, muQin, muSin);
    return eVal + pVal;
}

bool eos::update_s(double sin) { //update the t position (mu=0) based on input. Returns 1 if found, returns 0 if failed
    return update_s(sin, 0.0, 0.0, 0.0);
}

bool eos::update_s(double sin, double Bin, double Sin, double Qin) { //update the t and mu position based on input. Returns 1 if found, returns 0 if failed
    if (rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }
	///////////////////////////
	if (accept_nearest_neighbor)
		std::cout << "Should not have made it here!" << std::endl;
	if (discard_unsolvable_charge_densities)
		return false;//!!!!!!!!!!!!
	///////////////////////////
    double t0 = tbqsPosition[0];
    double mub0 = tbqsPosition[1];
    double muq0 = tbqsPosition[2];
    double mus0 = tbqsPosition[3];
    double t10 = t0*.2;
    double muB10 = mub0*.2;
    double muQ10 = muq0*.2;
    double muS10 = mus0*.2;

    //perturb T
    if(t0 + t10 > maxT) {
        tbqs(maxT - 1, mub0, muq0, mus0);
    } else {
        tbqs(t0 + t10, mub0, muq0, mus0);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }
    if(t0 - t10 < minT) {
        tbqs(minT + 1, mub0, muq0, mus0);
    } else {
        tbqs(t0 - t10, mub0, muq0, mus0);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }

    //perturb mub
    if(mub0 + muB10 > maxMuB) {
        tbqs(t0, maxMuB - 1, muq0, mus0);
    } else {
        tbqs(t0, mub0 + muB10, muq0, mus0);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }
    if(mub0 - muB10 < minMuB) {
        tbqs(t0, minMuB + 1, muq0, mus0);
    } else {
        tbqs(t0, mub0 - muB10, muq0, mus0);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }

    //perturn muq
    if(muq0 + muQ10 > maxMuQ) {
        tbqs(t0, mub0, maxMuQ - 1, mus0);
    } else {
        tbqs(t0, mub0, muq0 + muQ10, mus0);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }
    if(muq0 - muQ10 < minMuQ) {
        tbqs(t0, mub0, minMuQ + 1, mus0);
    } else {
        tbqs(t0, mub0, muq0 - muQ10, mus0);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }

    //perturb mus
    if(mus0 + muS10 > maxMuS) {
        tbqs(t0, mub0, muq0, maxMuS - 1);
    } else {
        tbqs(t0, mub0, muq0, mus0 + muS10);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }
    if(mus0 - muS10 < maxMuS) {
        tbqs(t0, mub0, muq0, minMuS + 1);
    } else {
        tbqs(t0, mub0, muq0, mus0 - muS10);
    }
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }

    //check mu = 0
    tbqs(t0, 0, 0, 0);
    if(rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }

    tbqs(t0, mub0, muq0, mus0);
    return false;
}



double eos::s_out(double ein) {   //update the t position (mu=0) based on input. Returns entropy if found, returns -1 if failed
    return s_out(ein, 0.0, 0.0, 0.0);
}

double eos::s_out(double ein, double Bin, double Sin, double Qin) {   //update the t and mu position based on input. Returns entropy if found, returns -1 if failed
    if (rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return entrVal;
    }
	///////////////////////////
	if (accept_nearest_neighbor)
		std::cout << "Should not have made it here!" << std::endl;
	if (discard_unsolvable_charge_densities)
		return -1.0;//!!!!!!!!!!!!
	///////////////////////////

    double t0 = tbqsPosition[0];
    double mub0 = tbqsPosition[1];
    double muq0 = tbqsPosition[2];
    double mus0 = tbqsPosition[3];
    double t10 = t0*.2;
    double muB10 = mub0*.2;
    double muQ10 = muq0*.2;
    double muS10 = mus0*.2;

    //perturb T
    if(t0 + t10 > maxT) {
        tbqs(maxT - 1, mub0, muq0, mus0);
    } else {
        tbqs(t0 + t10, mub0, muq0, mus0);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return entrVal;
    }
    if(t0 - t10 < minT) {
        tbqs(minT + 1, mub0, muq0, mus0);
    } else {
        tbqs(t0 - t10, mub0, muq0, mus0);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return entrVal;
    }

    //perturb mub
    if(mub0 + muB10 > maxMuB) {
        tbqs(t0, maxMuB - 1, muq0, mus0);
    } else {
        tbqs(t0, mub0 + muB10, muq0, mus0);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return entrVal;
    }
    if(mub0 - muB10 < minMuB) {
        tbqs(t0, minMuB + 1, muq0, mus0);
    } else {
        tbqs(t0, mub0 - muB10, muq0, mus0);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return entrVal;
    }

    //perturn muq
    if(muq0 + muQ10 > maxMuQ) {
        tbqs(t0, mub0, maxMuQ - 1, mus0);
    } else {
        tbqs(t0, mub0, muq0 + muQ10, mus0);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return entrVal;
    }
    if(muq0 - muQ10 < minMuQ) {
        tbqs(t0, mub0, minMuQ + 1, mus0);
    } else {
        tbqs(t0, mub0, muq0 - muQ10, mus0);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return entrVal;
    }

    //perturb mus
    if(mus0 + muS10 > maxMuS) {
        tbqs(t0, mub0, muq0, maxMuS - 1);
    } else {
        tbqs(t0, mub0, muq0, mus0 + muS10);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return entrVal;
    }
    if(mus0 - muS10 < maxMuS) {


        tbqs(t0, mub0, muq0, minMuS + 1);
    } else {
        tbqs(t0, mub0, muq0, mus0 - muS10);
    }
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return entrVal;
    }

    //check mu = 0
    tbqs(t0, 0, 0, 0);
    if(rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return entrVal;
    }

    tbqs(t0, mub0, muq0, mus0);
    return -1;
}


double eos::s_terms_T(double Tt) { //return entropy at a given temperature for muB = muS = muQ = 0
    tbqs(Tt, 0, 0, 0);
    return entrVal;
}


// UNCOMMENTED BY C. PLUMBERG
void eos::eosin(std::string type) {
}
double eos::A() {
    return w()-s()*dwds();
}


// confirm with Jaki
double eos::efreeze(double T_freeze_out_at_mu_eq_0) {
    tbqs(T_freeze_out_at_mu_eq_0, 0, 0, 0);
    return eVal;
}
double eos::sfreeze(double T_freeze_out_at_mu_eq_0) {
    return s_terms_T(T_freeze_out_at_mu_eq_0);
}




//struct to pass the target (E, rhoB, rhoQ, rhoS) into the rootfinder function
struct rootfinder_parameters {
    double eorEntGiven;          //these are the desired s and BSQ
    double rhoBGiven;
    double rhoQGiven;
    double rhoSGiven;
    rootfinder_parameters();
    rootfinder_parameters(
		double seteorEntGiven, double setRhoBGiven,
		double setRhoQGiven, double setRhoSGiven);
public:
    void set( double setEorEntGiven, double setRhoBGiven,
			  double setRhoQGiven, double setRhoSGiven);
};
//Default constructor to make the compiler happy. Should never be called
rootfinder_parameters::rootfinder_parameters() {}
//constructor which initializes all struct variables
rootfinder_parameters::rootfinder_parameters(
	double setEorEntGiven, double setRhoBGiven,
	double setRhoQGiven, double setRhoSGiven
	)
{
    eorEntGiven = setEorEntGiven;
    rhoBGiven = setRhoBGiven;
    rhoQGiven = setRhoQGiven;
    rhoSGiven = setRhoSGiven;
}
void rootfinder_parameters::set(
	double setEorEntGiven, double setRhoBGiven,
	double setRhoQGiven, double setRhoSGiven)
{
    eorEntGiven = setEorEntGiven;
    rhoBGiven = setRhoBGiven;
    rhoQGiven = setRhoQGiven;
    rhoSGiven = setRhoSGiven;

}

//helper function for the rootfinder. It provides the correct difference of s, rhoB, rhoQ, rhoS at a given (T, muB, muQ, muS) from the target
//used when rootfinder is given an entropy, baryon density, charge density, strangeness density
//x = (T, muB, muQ, muS), params = (sGiven, rhoBGiven, rhoQGiven, rhoSGiven), f becomes (s - sGiven, rhoB - rhoBGiven, rhoQ - rhoQGiven, rhoS - rhoSGiven)
int rootfinder_fsbqs(const gsl_vector *x, void *params, gsl_vector *f);
int rootfinder_fsbqs(const gsl_vector *x, void *params, gsl_vector *f) {
    //x contains the next (T, muB, muS) coordinate to test
    vector<double> tbqsToEval(4);
    tbqsToEval[0] = gsl_vector_get(x,0);
    tbqsToEval[1] = gsl_vector_get(x,1);	// convert x into densevector so it
    tbqsToEval[2] = gsl_vector_get(x,2);	// can be a BSpline evaluation point
    tbqsToEval[3] = gsl_vector_get(x,3);


    double entrGiven, rhoBGiven, rhoQGiven, rhoSGiven, entr, rhoB, rhoQ, rhoS;
    entrGiven = ((rootfinder_parameters*)params)->eorEntGiven;
    rhoBGiven = ((rootfinder_parameters*)params)->rhoBGiven;            //given variables contain the target point
    rhoQGiven = ((rootfinder_parameters*)params)->rhoQGiven;
    rhoSGiven = ((rootfinder_parameters*)params)->rhoSGiven;
{
	double phase_diagram_point[4] = {tbqsToEval[0]*197.3, tbqsToEval[1]*197.3,
					 tbqsToEval[3]*197.3, tbqsToEval[2]*197.3};	// NOTE: S <<-->> Q swapped!!!
	double densities_at_point[4];
	get_sBSQ_densities(phase_diagram_point, densities_at_point);
	entr = densities_at_point[0];
	rhoB = densities_at_point[1];
	rhoS = densities_at_point[2];
	rhoQ = densities_at_point[3];
	/*cout << "Check here: " 
		<< tbqsToEval(0)*197.3 << "   "
		<< tbqsToEval(1)*197.3 << "   "
		<< tbqsToEval(2)*197.3 << "   "
		<< tbqsToEval(3)*197.3 << "   "
		<< entr << "   " << entrGiven << "   "
		<< rhoB << "   " << rhoBGiven << "   "
		<< rhoS << "   " << rhoSGiven << "   "
		<< rhoQ << "   " << rhoQGiven << endl;*/

}

    gsl_vector_set(f, 0, (entr - entrGiven)); //f[0] contains (s(T,muB,muQ,muS) - sGiven)
    gsl_vector_set(f, 1, (rhoB - rhoBGiven)); //f[1] contains (rhoB(T,muB,muQ,muS) - rhoBGiven)
    gsl_vector_set(f, 2, (rhoQ - rhoQGiven)); //f[2] contains (rhoQ(T,muB,muQ,muS) - rhoQGiven)
    gsl_vector_set(f, 3, (rhoS - rhoSGiven)); //f[2] contains (rhoS(T,muB,muQ,muS) - rhoSGiven)

    return GSL_SUCCESS;
}

//helper function for the rootfinder. It provides the correct difference of E and rhoB at a given (T, muB, muQ, muS) from the target
//used when rootfinder is given an energy density and a baryon density
//x = (T, muB, muQ, muS), params = ((eGiven, rhoBGiven, rhoQGiven, rhoSGiven), f becomes (e - eGiven, rhoB - rhoBGiven, rhoQ - rhoQGiven, rhoS - rhoSGiven)
int rootfinder_febqs(const gsl_vector *x, void *params, gsl_vector *f);
int rootfinder_febqs(const gsl_vector *x, void *params, gsl_vector *f) {
    //x contains the next (T, muB, muQ, muS) coordinate to test
    vector<double> tbqsToEval(4);
    tbqsToEval[0] = gsl_vector_get(x,0);
    tbqsToEval[1] = gsl_vector_get(x,1);	// convert x into densevector so it
    tbqsToEval[2] = gsl_vector_get(x,2);	// can be a BSpline evaluation point
    tbqsToEval[3] = gsl_vector_get(x,3);

    double eGiven, rhoBGiven, rhoQGiven, rhoSGiven, e, rhoB, rhoQ, rhoS;
    eGiven = ((rootfinder_parameters*)params)->eorEntGiven;
    rhoBGiven = ((rootfinder_parameters*)params)->rhoBGiven;            //given variables contain the target point
    rhoQGiven = ((rootfinder_parameters*)params)->rhoQGiven;
    rhoSGiven = ((rootfinder_parameters*)params)->rhoSGiven;
{
	double phase_diagram_point[4] = {tbqsToEval[0]*197.3, tbqsToEval[1]*197.3,
					 tbqsToEval[3]*197.3, tbqsToEval[2]*197.3};	// NOTE: S <<-->> Q swapped!!!
	double densities_at_point[4];
	get_eBSQ_densities(phase_diagram_point, densities_at_point);
	e = densities_at_point[0]/197.3;
	rhoB = densities_at_point[1];
	rhoS = densities_at_point[2];
	rhoQ = densities_at_point[3];
	/*cout << "Check here: " 
		<< tbqsToEval(0)*197.3 << "   "
		<< tbqsToEval(1)*197.3 << "   "
		<< tbqsToEval(2)*197.3 << "   "
		<< tbqsToEval(3)*197.3 << "   "
		<< e*197.3 << "   " << eGiven*197.3 << "   "
		<< rhoB << "   " << rhoBGiven << "   "
		<< rhoS << "   " << rhoSGiven << "   "
		<< rhoQ << "   " << rhoQGiven << endl;
	*/
}

    gsl_vector_set(f, 0, (e - eGiven)); //f[0] contains (e(T,muB,muQ,muS) - eGiven)
    gsl_vector_set(f, 1, (rhoB - rhoBGiven)); //f[1] contains the (rhoB(T,muB,muQ,muS) - rhoBGiven)
    gsl_vector_set(f, 2, (rhoQ - rhoQGiven)); //f[2] contains the (rhoQ(T,muB,muQ,muS) - rhoQGiven)
    gsl_vector_set(f, 3, (rhoS - rhoSGiven)); //f[2] contains the (rho2(T,muB,muQ,muS) - rhoSGiven)

    return GSL_SUCCESS;
}


bool eos::rootfinder4D(double e_or_s_Given, int e_or_s_mode,
						double rhoBGiven, double rhoSGiven, double rhoQGiven,
						double error, size_t steps)
{
	if ( VERBOSE > 5 ) std::cout << __PRETTY_FUNCTION__ << e_or_s_Given << "   " << e_or_s_mode << "   " << rhoBGiven << "   " << rhoSGiven << "   " << rhoQGiven << "   " << error << "   " << steps << std::endl;

    //declare x = (T, muB, muQ, muS)
    gsl_vector *x = gsl_vector_alloc(4);

	// use NMN method to estimate where to start the rootfinder
	// ( returns estimates in units of MeV )
	vector<double> T_muB_muQ_muS_estimates;
	constexpr bool use_normalized_trees = true;
	if ( e_or_s_mode==1 )
		e_delaunay.get_NMN_coordinates(
					{e_or_s_Given*197.3, rhoBGiven, rhoSGiven, rhoQGiven},
					T_muB_muQ_muS_estimates, use_normalized_trees );
	else
		entr_delaunay.get_NMN_coordinates(
					{e_or_s_Given, rhoBGiven, rhoSGiven, rhoQGiven},
					T_muB_muQ_muS_estimates, use_normalized_trees );

	// set GSL vector with best initial guess we can
	/*std::cout << "Closest neighbor found to be:";
	for (int iCoord = 0; iCoord < 4; iCoord++)
	{
		std::cout << "   " << T_muB_muQ_muS_estimates[iCoord];
		gsl_vector_set(x, iCoord, T_muB_muQ_muS_estimates[iCoord]/197.3);
	}
	std::cout << std::endl;*/

    gsl_vector_set(x, 0, T());
    gsl_vector_set(x, 1, muB());
    gsl_vector_set(x, 2, muQ());
    gsl_vector_set(x, 3, muS());
    /*gsl_vector_set(x, 0, 500.0/197.3);
    gsl_vector_set(x, 1, 0.0/197.3);
    gsl_vector_set(x, 2, 0.0/197.3);
    gsl_vector_set(x, 3, 0.0/197.3);*/

    //initialize the rootfinder equation to the correct variable quantities
    bool isEntropy = false;
    if(e_or_s_mode == 0) {
        isEntropy = true;
    }
    rootfinder_parameters p;
	if(isEntropy) {
		p.set( e_or_s_Given, rhoBGiven, rhoQGiven, rhoSGiven);
	} else {
		p.set( e_or_s_Given, rhoBGiven, rhoQGiven, rhoSGiven);
	}

    //initialize multiroot solver
    gsl_multiroot_fsolver *solver;
    gsl_multiroot_function f;
    if(isEntropy) {
        f.f = &rootfinder_fsbqs;
    } else {
        f.f = &rootfinder_febqs;
    }
    f.n = 4;
    f.params = &p;

	/*if ( e_or_s_mode == 1 and VERBOSE > 5 )
		std::cout << std::endl
			<< "=============================================="
			<< std::endl << "Input (e,B,Q,S): "
			<< e_or_s_Given*0.1973 << "   "
			<< rhoBGiven << "   "
			<< rhoQGiven << "   "
			<< rhoSGiven << std::endl;*/

    solver = gsl_multiroot_fsolver_alloc(TYPE, 4);
    gsl_multiroot_fsolver_set(solver, &f, x);

    int status;
    size_t iter = 0;
	double previous_solver_step[4];

    do
	{
		for (int iSolverElem = 0; iSolverElem < 4; iSolverElem++)
			previous_solver_step[iSolverElem] = gsl_vector_get(solver->x, iSolverElem);

        ++iter;
        status = gsl_multiroot_fsolver_iterate(solver);

        if(VERBOSE > 5 && status)
		{
			if ( status == GSL_EBADFUNC && VERBOSE > 5 )
				std::cout << "Error: something went to +/-Inf or NaN!" << std::endl;
			else if ( status == GSL_ENOPROG && VERBOSE > 5 )
				std::cout << "Error: not making enough progress!" << std::endl;
            //break if the rootfinder gets stuck
			break;
        }

		//break if the rootfinder goes out of bounds
        if(gsl_vector_get(solver->x, 0) < minT)
		{
			if ( VERBOSE > 5 )
				std::cout << "Error: out-of-bounds (T < minT)!" << std::endl;
			status = -10;
			break;
        }
		else if(gsl_vector_get(solver->x, 0) > maxT)
		{
			if ( VERBOSE > 5 )
				std::cout << "Error: out-of-bounds (T > maxT)!" << std::endl;
			status = -10;
			break;
        }
		else if (gsl_vector_get(solver->x, 1) < minMuB)
		{
			if ( VERBOSE > 5 )
				std::cout << "Error: out-of-bounds (MuB < minMuB)!" << std::endl;
			status = -10;
			break;
        }
		else if (gsl_vector_get(solver->x, 1) > maxMuB)
		{
			if ( VERBOSE > 5 )
				std::cout << "Error: out-of-bounds (MuB > maxMuB)!" << std::endl;
			status = -10;
			break;
        }
		else if (gsl_vector_get(solver->x, 2) < minMuQ)
		{
			if ( VERBOSE > 5 )
				std::cout << "Error: out-of-bounds (MuQ < minMuQ)!" << std::endl;
			status = -10;
			break;
        }
		else if (gsl_vector_get(solver->x, 2) > maxMuQ)
		{
			if ( VERBOSE > 5 )
				std::cout << "Error: out-of-bounds (MuQ > maxMuQ)!" << std::endl;
			status = -10;
			break;
        }
		else if (gsl_vector_get(solver->x, 3) < minMuS)
		{
			if ( VERBOSE > 5 )
				std::cout << "Error: out-of-bounds (MuS < minMuS)!" << std::endl;
			status = -10;
			break;
        }
		else if (gsl_vector_get(solver->x, 3) > maxMuS)
		{
			if ( VERBOSE > 5 )
				std::cout << "Error: out-of-bounds (MuS > maxMuS)!" << std::endl;
			status = -10;
			break;
        }

        status = gsl_multiroot_test_residual(solver->f, error);

    } while (status == GSL_CONTINUE && iter < steps);

//std::cout << "Exited GSL loop" << endl;

    bool found = true; //to return variable
    if ( iter >= steps || status != 0 )
	{
		if ( status == GSL_EBADFUNC )
			std::cout << "Error: something went to +/-Inf or NaN!" << std::endl;
		else if ( status == GSL_ENOPROG )
			std::cout << "Error: not making enough progress!" << std::endl;
		else
			std::cout << "Check: " << iter << "   " << steps << "   " << status << std::endl;
        found = false;
    }

//std::cout << "found = " << found << endl;


    if ( found )
	{
        tbqs( gsl_vector_get(solver->x, 0),
			  gsl_vector_get(solver->x, 1),
			  gsl_vector_get(solver->x, 2),
			  gsl_vector_get(solver->x, 3) );    //set T, muB, muQ, muS
		/*cout << "GSL found solution!";
		for (int iSol = 0; iSol <4; iSol++)
			cout << "   " << gsl_vector_get(solver->x, iSol);
		cout << endl;*/
    }
	else if ( accept_nearest_neighbor )	// just return nearest neighbor instead
	{
		std::cout << "Entered to nearest neighbor tests" << std::endl;
		std::cout << "found = " << found << std::endl;

		// use unnormalized distances to estimate neighbor (reset from above)
		if ( e_or_s_mode==1 )
			e_delaunay.get_NMN_coordinates(
						{e_or_s_Given*197.3, rhoBGiven, rhoSGiven, rhoQGiven},
						T_muB_muQ_muS_estimates, false );
		else
			entr_delaunay.get_NMN_coordinates(
						{e_or_s_Given, rhoBGiven, rhoSGiven, rhoQGiven},
						T_muB_muQ_muS_estimates, false );

		// set final solver location
		int which_neighbor_closest = -1;

		double Tfinal = 0.0, muBfinal = 0.0, muQfinal = 0.0, muSfinal = 0.0;
		if (iter >= steps && status == 0)	// no reported problems, just ran out of steps
		{
			Tfinal   = gsl_vector_get(solver->x, 0)*197.3;
			muBfinal = gsl_vector_get(solver->x, 1)*197.3;
			muQfinal = gsl_vector_get(solver->x, 2)*197.3;
			muSfinal = gsl_vector_get(solver->x, 3)*197.3;
		}
		else
		{
			Tfinal   = previous_solver_step[0]*197.3;
			muBfinal = previous_solver_step[1]*197.3;
			muQfinal = previous_solver_step[2]*197.3;
			muSfinal = previous_solver_step[3]*197.3;
		}

		double inputDensities[4] = {e_or_s_Given, rhoBGiven, rhoSGiven, rhoQGiven};
		double finalDensities[4], neighbor_estimate_densities[4];
		// swap muQ <--> muS!!!
		double final_phase_diagram_point[4] = {Tfinal, muBfinal, muSfinal, muQfinal};
		// swap muQ <--> muS!!!
		double neighbor_estimate_point[4]
				= {T_muB_muQ_muS_estimates[0], T_muB_muQ_muS_estimates[1],
					T_muB_muQ_muS_estimates[3], T_muB_muQ_muS_estimates[2]};
		std::cout << "final_phase_diagram_point =";
		for (int iii = 0; iii < 4; iii++)
			std::cout << "   " << final_phase_diagram_point[iii];
		std::cout << std::endl;
		std::cout << "neighbor_estimate_point =";
		for (int iii = 0; iii < 4; iii++)
			std::cout << "   " << neighbor_estimate_point[iii];
		std::cout << std::endl;

		if ( isEntropy )
		{
			get_sBSQ_densities(final_phase_diagram_point, finalDensities);
			get_sBSQ_densities(neighbor_estimate_point, neighbor_estimate_densities);
			which_neighbor_closest
				= static_cast<int>(
					entr_delaunay.unnormalized_d2( inputDensities, neighbor_estimate_densities )
					< entr_delaunay.unnormalized_d2( inputDensities, finalDensities ) );
			std::cout << "Check separations: "
					<< entr_delaunay.unnormalized_d2( inputDensities,
													neighbor_estimate_densities ) << "   "
					<< entr_delaunay.unnormalized_d2( inputDensities, finalDensities ) << "   "
					<< which_neighbor_closest << std::endl;
		}
		else
		{
			get_eBSQ_densities(final_phase_diagram_point, finalDensities);
			get_eBSQ_densities(neighbor_estimate_point, neighbor_estimate_densities);
			which_neighbor_closest
				= static_cast<int>(
					e_delaunay.unnormalized_d2( inputDensities, neighbor_estimate_densities )
					< e_delaunay.unnormalized_d2( inputDensities, finalDensities ) );
			std::cout << "Check separations: "
					<< e_delaunay.unnormalized_d2( inputDensities,
													neighbor_estimate_densities ) << "   "
					<< e_delaunay.unnormalized_d2( inputDensities, finalDensities ) << "   "
					<< which_neighbor_closest << std::endl;
		}

		std::cout << "Made it to line = " << __LINE__ << std::endl;

		// set (T, muB, muQ, muS) based on which point is closest to input point
		// (SWAP S <--> Q AGAIN)
		if ( which_neighbor_closest == 0 )
		{
			tbqs( final_phase_diagram_point[0]/197.3, final_phase_diagram_point[1]/197.3,
				  final_phase_diagram_point[3]/197.3, final_phase_diagram_point[2]/197.3 );
		cout << "final_phase_diagram_point: ";
		for (int iSol = 0; iSol <4; iSol++)
			cout << "   " << final_phase_diagram_point[iSol] / 197.3;
		cout << endl;

		}
		else if ( which_neighbor_closest == 1 )
		{
			tbqs( neighbor_estimate_point[0]/197.3, neighbor_estimate_point[1]/197.3,
				  neighbor_estimate_point[3]/197.3, neighbor_estimate_point[2]/197.3 );
		cout << "neighbor_estimate_point!";
		for (int iSol = 0; iSol < 4; iSol++)
			cout << "   " << neighbor_estimate_point[iSol] / 197.3;
		cout << endl;

		}
		else
		{
			std::cerr << "Bad value: which_neighbor_closest = "
					<< which_neighbor_closest << std::endl;
			exit(-8);
		}

		std::cout << "Setting found --> true" << std::endl;
		found = true;
	}
	/*else if (!discard_unsolvable_charge_densities)
	{
		std::cerr << "ERROR: you still need a back-up plan!" << std::endl;
		exit(-8);
	}*/


    //memory deallocation
    gsl_multiroot_fsolver_free(solver);
    gsl_vector_free(x);
    return found;
}


