#include "eos.h"

#include "../splinter/include/datatable.h"
#include "../splinter/include/bspline.h"
#include "../splinter/include/bsplinebuilder.h"
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

using namespace SPLINTER;

// Compile:         gcc eos4D.cpp -c -I /usr/include/eigen3 -Lsplinter/build -lm -lgsl -lgslcblas -lstdc++ -lsplinter-3-0


//EoS constructor. Builds the splines of degree "degree" for each quantitiy and initializes the position at (30,0,0,0)
eos::eos(string quantityFile, string derivFile, int degree) : pSpline(4), entrSpline(4), bSpline(4), sSpline(4), qSpline(4), eSpline(4), cs2Spline(4), db2Spline(4), dq2Spline(4), ds2Spline(4), dt2Spline(4), dbdqSpline(4), dbdsSpline(4), dtdbSpline(4), dqdsSpline(4), dtdqSpline(4), dtdsSpline(4), tbqsPosition(4) {
    init(quantityFile, derivFile, degree);
}

//EoS default constructor. This function exists to satisfy the compiler
//This function should never be called unless init is called directly afterward
eos::eos() : pSpline(4), entrSpline(4), bSpline(4), sSpline(4), qSpline(4), eSpline(4), cs2Spline(4), db2Spline(4), dq2Spline(4), ds2Spline(4), dt2Spline(4), dbdqSpline(4), dbdsSpline(4), dtdbSpline(4), dqdsSpline(4), dtdqSpline(4), dtdsSpline(4), tbqsPosition(4) {}

void eos::init(string quantityFile, string derivFile, int degree)
{
	VERBOSE = 5;

	cout << "Initializing EoS C library" << endl;
	initialize("/projects/jnorhos/BSQ/EoS_BQS_Derivatives/Coefficients_Parameters.dat");

	std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
	init_with_txt(quantityFile, derivFile, degree);

	cout << "Initialize Delaunay interpolators" << endl;
	e_delaunay.init(quantityFile, 0);		// 0 - energy density
	entr_delaunay.init(quantityFile, 1);	// 1 - entropy density

	return;
}

void eos::init_with_txt(string quantityFile, string derivFile, int degree)
{
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    std::ifstream dataFile;
    std::ifstream derFile;
    dataFile.open(quantityFile);
    derFile.open(derivFile);

    DataTable psamples, entrsamples, bsamples, ssamples, qsamples, esamples, cs2samples;
    DataTable db2samples, ds2samples, dq2samples, dt2samples, dbdssamples, dbdqsamples, dqdssamples, dtdssamples, dtdqsamples, dtdbsamples;

    double tit, muBit, muQit, muSit, pit, entrit, bit, sit, qit, eit, cs2it;
    double db2it, dq2it, ds2it, dt2it, dbdqit, dbdsit, dqdsit, dtdbit, dtdsit, dtdqit;
    vector<double> toAdd;

    int count = 0;
    double hbarc = 197.327;
    while (dataFile >> tit >> muBit >> muQit >> muSit >> pit >> entrit >> bit >> sit >> qit >> eit >> cs2it) {
        derFile >> tit >> muBit >> muQit >> muSit >> db2it >> dq2it >> ds2it >> dbdqit >> dbdsit >> dqdsit >> dtdbit >> dtdqit >> dtdsit >> dt2it;  //read data from files

		// Christopher Plumberg:
		// put T and mu_i in units of 1/fm
		tit   /= hbarc;
		muBit /= hbarc;
		muSit /= hbarc;
		muQit /= hbarc;

        if(count++ == 0) {
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
        if(maxT < tit) {
            maxT = tit;
        }
        if(minT > tit) {
            minT = tit;
        }
        if(maxMuB < muBit) {
            maxMuB = muBit;
        }
        if(minMuB > muBit) {
            minMuB = muBit;
        }
        if(maxMuQ < muQit) {
            maxMuQ = muQit;
        }
        if(minMuQ > muQit) {
            minMuQ = muQit;
        }
        if(maxMuS < muSit) {
            maxMuS = muSit;
        }
        if(minMuS > muSit) {
            minMuS = muSit;
        }

        toAdd.push_back(tit);
        toAdd.push_back(muBit);
        toAdd.push_back(muQit);
        toAdd.push_back(muSit);

		/*
        pit = pit*(tit*tit*tit*tit)/(hbarc*hbarc*hbarc);
        entrit = entrit*(tit*tit*tit)/(hbarc*hbarc*hbarc);
        bit = bit*(tit*tit*tit)/(hbarc*hbarc*hbarc);
        sit = sit*(tit*tit*tit)/(hbarc*hbarc*hbarc);        //!!!!!convert to MeV and fm units
        qit = qit*(tit*tit*tit)/(hbarc*hbarc*hbarc);		//!!!!! --> NOW fm only for use in hydro!
        eit = eit*(tit*tit*tit*tit)/(hbarc*hbarc*hbarc);
		*/

		// USE FM IN HYDRO
        pit = pit*(tit*tit*tit*tit);
        entrit = entrit*(tit*tit*tit);
        bit = bit*(tit*tit*tit);
        sit = sit*(tit*tit*tit);
        qit = qit*(tit*tit*tit);
        eit = eit*(tit*tit*tit*tit);
		

        psamples.addSample(toAdd, pit);
        entrsamples.addSample(toAdd, entrit);
        bsamples.addSample(toAdd, bit);
        ssamples.addSample(toAdd, sit);
        qsamples.addSample(toAdd, qit);
        esamples.addSample(toAdd, eit);
        cs2samples.addSample(toAdd, cs2it);
        db2samples.addSample(toAdd, db2it);
        dq2samples.addSample(toAdd, dq2it);     //add datapoint to table for spline builder
        ds2samples.addSample(toAdd, ds2it);
        dbdqsamples.addSample(toAdd, dbdqit);
        dbdssamples.addSample(toAdd, dbdsit);
        dqdssamples.addSample(toAdd, dqdsit);
        dtdbsamples.addSample(toAdd, dtdbit);
        dtdqsamples.addSample(toAdd, dtdqit);
        dtdssamples.addSample(toAdd, dtdsit);
        dt2samples.addSample(toAdd, dt2it);
        toAdd.clear();
    }

    dataFile.close();
    derFile.close();

	std::cout << "Finished reading in thermodynamic data files!" << std::endl;

	std::cout << "Building pspline..." << std::endl;
    pSpline = BSpline::Builder(psamples).degree(degree).build();
	std::cout << "Building entrSpline..." << std::endl;
    entrSpline = BSpline::Builder(entrsamples).degree(degree).build();
	std::cout << "Building bSpline..." << std::endl;
    bSpline = BSpline::Builder(bsamples).degree(degree).build();
	std::cout << "Building sSpline..." << std::endl;
    sSpline = BSpline::Builder(ssamples).degree(degree).build();
	std::cout << "Building qSpline..." << std::endl;
    qSpline = BSpline::Builder(qsamples).degree(degree).build();
	std::cout << "Building eSpline..." << std::endl;
    eSpline = BSpline::Builder(esamples).degree(degree).build();
	std::cout << "Building cs2Spline..." << std::endl;
    cs2Spline = BSpline::Builder(cs2samples).degree(degree).build();
	std::cout << "Building db2Spline..." << std::endl;
    db2Spline = BSpline::Builder(db2samples).degree(degree).build();
	std::cout << "Building dq2Spline..." << std::endl;
    dq2Spline = BSpline::Builder(dq2samples).degree(degree).build();
 	std::cout << "Building ds2Spline..." << std::endl;
    ds2Spline = BSpline::Builder(ds2samples).degree(degree).build();        //make splines from table
	std::cout << "Building dbdqSpline..." << std::endl;
    dbdqSpline = BSpline::Builder(dbdqsamples).degree(degree).build();
	std::cout << "Building dbdsSpline..." << std::endl;
    dbdsSpline = BSpline::Builder(dbdssamples).degree(degree).build();
	std::cout << "Building dqdsSpline..." << std::endl;
    dqdsSpline = BSpline::Builder(dqdssamples).degree(degree).build();
	std::cout << "Building dtdbSpline..." << std::endl;
    dtdbSpline = BSpline::Builder(dtdbsamples).degree(degree).build();
	std::cout << "Building dtdqSpline..." << std::endl;
    dtdqSpline = BSpline::Builder(dtdqsamples).degree(degree).build();
	std::cout << "Building dtdsSpline..." << std::endl;
    dtdsSpline = BSpline::Builder(dtdssamples).degree(degree).build();
	std::cout << "Building dt2Spline..." << std::endl;
    dt2Spline = BSpline::Builder(dt2samples).degree(degree).build();

	// initialize tbqsPosition to something...
	std::cout << "Initializing tbqsPosition...\n";
	for (int iTBQS = 0; iTBQS < 4; iTBQS++) tbqsPosition(iTBQS) = 1.0;

	std::cout << "Check TBQS: ";
	for (int iTBQS = 0; iTBQS < 4; iTBQS++) std::cout << tbqsPosition(iTBQS) << "   ";	
	std::cout << std::endl;

	//std::cout << "Check alternate: "
	//			<< T() << "   " << muB() << "   "
	//			<< muQ() << "   " << muS() << std::endl;

	std::cout << "All initializations finished!" << std::endl;

    return;
}

void eos::tbqs(double setT, double setmuB, double setmuQ, double setmuS) {
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;

    if(setT < minT || setT > maxT) {
        std::cout << "T = " << setT << " is out of range. Valid values are between [" << minT << "," << maxT << "]" << std::endl;
        return;
    }
    if(setmuB < minMuB || setmuB > maxMuB) {
        std::cout << "muB = " << setmuB << " is out of range. Valid values are between [" << minMuB << "," << maxMuB << "]" << std::endl;
        return;
    }
    if(setmuQ < minMuQ || setmuQ > maxMuQ) {
        std::cout << "muQ = " << setmuQ << " is out of range. Valid values are between [" << minMuQ << "," << maxMuQ << "]" << std::endl;
        return;
    }
    if(setmuS < minMuS || setmuS > maxMuS) {
        std::cout << "muS = " << setmuS << " is out of range. Valid values are between [" << minMuS << "," << maxMuS << "]" << std::endl;
        return;
    }
    tbqsPosition(0) = setT;
    tbqsPosition(1) = setmuB;
    tbqsPosition(2) = setmuQ;
    tbqsPosition(3) = setmuS;

    pVal = pSpline.eval(tbqsPosition);
    BVal = bSpline.eval(tbqsPosition);
    SVal = sSpline.eval(tbqsPosition);
    QVal = qSpline.eval(tbqsPosition);
    eVal = eSpline.eval(tbqsPosition);
    cs2Val = cs2Spline.eval(tbqsPosition);
    db2 = db2Spline.eval(tbqsPosition);
    ds2 = ds2Spline.eval(tbqsPosition);
    dq2 = dq2Spline.eval(tbqsPosition);
    dt2 = dt2Spline.eval(tbqsPosition);
    dbdq = dbdqSpline.eval(tbqsPosition);
    dbds = dbdsSpline.eval(tbqsPosition);
    dsdq = dqdsSpline.eval(tbqsPosition);
    dtdb = dtdbSpline.eval(tbqsPosition);
    dtds = dtdsSpline.eval(tbqsPosition);
    dtdq = dtdqSpline.eval(tbqsPosition);

    entrVal = (eVal + pVal - setmuB*BVal - setmuQ*QVal - setmuS*SVal)/setT;
}


double eos::T() {
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << ": " << tbqsPosition(0) << std::endl;
    return tbqsPosition(0);
}

double eos::muB() {
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << ": " << tbqsPosition(1) << std::endl;
    return tbqsPosition(1);
}

double eos::muQ() {
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << ": " << tbqsPosition(2) << std::endl;
    return tbqsPosition(2);
}

double eos::muS() {
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << ": " << tbqsPosition(3) << std::endl;
    return tbqsPosition(3);
}




double eos::p() {
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return pVal;
}

double eos::s() {
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return entrVal;
}

double eos::B() {
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return BVal;
}

double eos::S() {
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return SVal;
}

double eos::Q() {
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return QVal;
}

double eos::e() {
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return eVal;
}

double eos::cs2() {
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return cs2Val;
}


double eos::w() {
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return eVal + pVal;
}

double eos::dwds()
{
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;

	double charge_terms	/*if charge densities are not all zero*/
			= ( abs(BVal)>1e-10 || abs(SVal)>1e-10 || abs(QVal)>1e-10 ) ?
			  BVal/dentr_dmub() + QVal/dentr_dmuq() + SVal/dentr_dmus() : 0.0;

    return T() + entrVal/dentr_dt() + charge_terms;
}

double eos::dwdB()
{
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;

	double charge_terms	/*if charge densities are not all zero*/
			= ( abs(BVal)>1e-10 || abs(SVal)>1e-10 || abs(QVal)>1e-10 ) ?
			  entrVal/db_dt() + BVal/db_dmub() + QVal/db_dmuq() + SVal/db_dmus() : 0.0;

    return muB() + charge_terms;
}

double eos::dwdS()
{
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;

	double charge_terms	/*if charge densities are not all zero*/
			= ( abs(BVal)>1e-10 || abs(SVal)>1e-10 || abs(QVal)>1e-10 ) ?
			  entrVal/ds_dt() + BVal/ds_dmub() + QVal/ds_dmuq() + SVal/ds_dmus() : 0.0;

    return muS() + charge_terms;
}

double eos::dwdQ()
{
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;

	double charge_terms	/*if charge densities are not all zero*/
			= ( abs(BVal)>1e-10 || abs(SVal)>1e-10 || abs(QVal)>1e-10 ) ?
			  entrVal/dq_dt() + BVal/dq_dmub() + QVal/dq_dmuq() + SVal/dq_dmus() : 0.0;

    return muQ() + charge_terms;
}

double eos::dentr_dt() {
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return calc_term_1();
}

double eos::dentr_dmub() {
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return calc_term_2("b");
}

double eos::dentr_dmuq() {
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return calc_term_2("q");
}

double eos::dentr_dmus() {
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return calc_term_2("s");
}

double eos::db_dt() {
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return calc_term_3("b");
}

double eos::db_dmub() {
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return calc_term_4("b","b");
}

double eos::db_dmuq() {
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return calc_term_4("b","q");
}

double eos::db_dmus() {
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return calc_term_4("b","s");
}

double eos::ds_dt() {
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return calc_term_3("s");
}

double eos::ds_dmub() {
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return calc_term_4("s","b");
}

double eos::ds_dmuq() {
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return calc_term_4("s","q");
}

double eos::ds_dmus() {
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return calc_term_4("s","s");
}

double eos::dq_dt() {
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return calc_term_3("q");
}

double eos::dq_dmub() {
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return calc_term_4("q","b");
}

double eos::dq_dmuq() {
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return calc_term_4("q","q");
}

double eos::dq_dmus() {
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return calc_term_4("q","s");
}

double eos::calc_term_1() {
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    gsl_vector *v = gsl_vector_alloc(3);
    gsl_matrix *m = gsl_matrix_alloc(3,3);

    gsl_vector_set(v,0,dtdb);
    gsl_vector_set(v,1,dtds);
    gsl_vector_set(v,2,dtdq);

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
	if ( VERBOSE > 1 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
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
				<< 197.327*T() << "   " << 197.327*muB() << "   "
				<< 197.327*muS() << "   " << 197.327*muQ() << endl << endl;

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

		cout << "m=" << endl;
		for (int ii = 0; ii < 3; ii++)
		{
			for (int jj = 0; jj < 3; jj++)
				cout << gsl_matrix_get(m, ii, jj) << "   ";
			cout << endl;
		}
		cout << "minv=" << endl;
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
    gsl_blas_dgemv(CblasNoTrans,1,m,b,0,y);

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
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    Aout=w()-entrVal*dwds();

    return Aout;
}

double eos::cs2out(double Tt) {  //return cs2 given t and mu's=0
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    tbqs(Tt, 0.0, 0.0, 0.0);
    return cs2Val;
}

double eos::cs2out(double Tt, double muBin, double muQin, double muSin) {  //return cs2 given t and mu's
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    tbqs(Tt, muBin, muQin, muSin);
    return cs2Val;
}

double eos::wfz(double Tt) {   // return e + p for tbqs
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    tbqs(Tt, 0.0, 0.0, 0.0);
    return eVal + pVal;
}

double eos::wfz(double Tt, double muBin, double muQin, double muSin) {   // return e + p for tbqs
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    tbqs(Tt, muBin, muQin, muSin);
    return eVal + pVal;
}

bool eos::update_s(double sin) { //update the t position (mu=0) based on input. Returns 1 if found, returns 0 if failed
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return update_s(sin, 0.0, 0.0, 0.0);
}

bool eos::update_s(double sin, double Bin, double Sin, double Qin) { //update the t and mu position based on input. Returns 1 if found, returns 0 if failed
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    if (rootfinder4D(sin, 0, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return true;
    }
    double t0 = tbqsPosition(0);
    double mub0 = tbqsPosition(1);
    double muq0 = tbqsPosition(2);
    double mus0 = tbqsPosition(3);
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
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return s_out(ein, 0.0, 0.0, 0.0);
}

double eos::s_out(double ein, double Bin, double Sin, double Qin) {   //update the t and mu position based on input. Returns entropy if found, returns -1 if failed
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    if (rootfinder4D(ein, 1, Bin, Sin, Qin, TOLERANCE, STEPS)) {
        return entrVal;
    }

    double t0 = tbqsPosition(0);
    double mub0 = tbqsPosition(1);
    double muq0 = tbqsPosition(2);
    double mus0 = tbqsPosition(3);
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
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    tbqs(Tt, 0, 0, 0);
    return entrVal;
}


// UNCOMMENTED BY C. PLUMBERG
void eos::eosin(std::string type) {
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
}
double eos::A() {
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return 0;
}


// confirm with Jaki
double eos::efreeze(double T_freeze_out_at_mu_eq_0) {
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    tbqs(T_freeze_out_at_mu_eq_0, 0, 0, 0);
    return eVal;
}
double eos::sfreeze(double T_freeze_out_at_mu_eq_0) {
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return s_terms_T(T_freeze_out_at_mu_eq_0);
}




//struct to pass the target (E, rhoB, rhoQ, rhoS) into the rootfinder function
struct rootfinder_parameters {
    double eorEntGiven;          //these are the desired s and BSQ
    double rhoBGiven;
    double rhoQGiven;
    double rhoSGiven;
    BSpline eorEntSpline;        //the splines that contain interpolations over s, BSQ
    BSpline rhoBSpline;
    BSpline rhoQSpline;
    BSpline rhoSSpline;
    rootfinder_parameters();
    rootfinder_parameters(double seteorEntGiven, double setRhoBGiven, double setRhoQGiven, double setRhoSGiven, BSpline setEntrSpline, BSpline setRhoBSPLine, BSpline setRhoQSpline, BSpline setRhoSSpline);
public:
    void set(double setEorEntGiven, double setRhoBGiven, double setRhoQGiven, double setRhoSGiven, BSpline setEntrSpline, BSpline setRhoBSpline, BSpline setRhoQSpline, BSpline setRhoSSpline);
};
//Default constructor to make the compiler happy. Should never be called
rootfinder_parameters::rootfinder_parameters() : eorEntSpline(4), rhoBSpline(4), rhoQSpline(4), rhoSSpline(4) {}
//constructor which initializes all struct variables
rootfinder_parameters::rootfinder_parameters(double setEorEntGiven, double setRhoBGiven, double setRhoQGiven, double setRhoSGiven, BSpline setEorEntSpline, BSpline setRhoBSpline, BSpline setRhoQSpline, BSpline setRhoSSpline) : eorEntSpline(4), rhoBSpline(4), rhoQSpline(4), rhoSSpline(4)
{
    eorEntGiven = setEorEntGiven;
    rhoBGiven = setRhoBGiven;
    rhoQGiven = setRhoQGiven;
    rhoSGiven = setRhoSGiven;
    eorEntSpline = setEorEntSpline;
    rhoBSpline = setRhoBSpline;
    rhoQSpline = setRhoQSpline;
    rhoSSpline = setRhoSSpline;
}
void rootfinder_parameters::set(double setEorEntGiven, double setRhoBGiven, double setRhoQGiven, double setRhoSGiven, BSpline setEorEntSpline, BSpline setRhoBSpline, BSpline setRhoQSpline, BSpline setRhoSSpline) {
    eorEntGiven = setEorEntGiven;
    rhoBGiven = setRhoBGiven;
    rhoQGiven = setRhoQGiven;
    rhoSGiven = setRhoSGiven;
    eorEntSpline = setEorEntSpline;
    rhoBSpline = setRhoBSpline;
    rhoQSpline = setRhoQSpline;
    rhoSSpline = setRhoSSpline;
}

//helper function for the rootfinder. It provides the correct difference of s, rhoB, rhoQ, rhoS at a given (T, muB, muQ, muS) from the target
//used when rootfinder is given an entropy, baryon density, charge density, strangeness density
//x = (T, muB, muQ, muS), params = (sGiven, rhoBGiven, rhoQGiven, rhoSGiven), f becomes (s - sGiven, rhoB - rhoBGiven, rhoQ - rhoQGiven, rhoS - rhoSGiven)
int rootfinder_fsbqs(const gsl_vector *x, void *params, gsl_vector *f);
int rootfinder_fsbqs(const gsl_vector *x, void *params, gsl_vector *f) {
    //x contains the next (T, muB, muS) coordinate to test
    DenseVector tbqsToEval(4);
    tbqsToEval(0) = gsl_vector_get(x,0);
    tbqsToEval(1) = gsl_vector_get(x,1);      //convert x into densevector so it can be a BSpline evaluation point
    tbqsToEval(2) = gsl_vector_get(x,2);
    tbqsToEval(3) = gsl_vector_get(x,3);


    double entrGiven, rhoBGiven, rhoQGiven, rhoSGiven, entr, rhoB, rhoQ, rhoS;
    entrGiven = ((rootfinder_parameters*)params)->eorEntGiven;
    rhoBGiven = ((rootfinder_parameters*)params)->rhoBGiven;            //given variables contain the target point
    rhoQGiven = ((rootfinder_parameters*)params)->rhoQGiven;
    rhoSGiven = ((rootfinder_parameters*)params)->rhoSGiven;
    entr = (((rootfinder_parameters*)params)->eorEntSpline).eval(tbqsToEval);    //s, rhoB, rhoQ, rhoS contain the current point
    rhoB = (((rootfinder_parameters*)params)->rhoBSpline).eval(tbqsToEval);
    rhoQ = (((rootfinder_parameters*)params)->rhoQSpline).eval(tbqsToEval);
    rhoS = (((rootfinder_parameters*)params)->rhoSSpline).eval(tbqsToEval);

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
    DenseVector tbqsToEval(4);
    tbqsToEval(0) = gsl_vector_get(x,0);
    tbqsToEval(1) = gsl_vector_get(x,1);      //convert x into densevector so it can be a BSpline evaluation point
    tbqsToEval(2) = gsl_vector_get(x,2);
    tbqsToEval(3) = gsl_vector_get(x,3);

    double eGiven, rhoBGiven, rhoQGiven, rhoSGiven, e, rhoB, rhoQ, rhoS;
    eGiven = ((rootfinder_parameters*)params)->eorEntGiven;
    rhoBGiven = ((rootfinder_parameters*)params)->rhoBGiven;            //given variables contain the target point
    rhoQGiven = ((rootfinder_parameters*)params)->rhoQGiven;
    rhoSGiven = ((rootfinder_parameters*)params)->rhoSGiven;
if (true)
{
	double phase_diagram_point[4] = {tbqsToEval(0)*197.327, tbqsToEval(1)*197.327,
					 tbqsToEval(3)*197.327, tbqsToEval(2)*197.327};	// NOTE: S <<-->> Q swapped!!!
	double densities_at_point[4];
	get_densities(phase_diagram_point, densities_at_point);
	e = densities_at_point[0]/197.327;
	rhoB = densities_at_point[1];
	rhoS = densities_at_point[2];
	rhoQ = densities_at_point[3];
	/*cout << "Check here: " 
		<< tbqsToEval(0)*197.327 << "   "
		<< tbqsToEval(1)*197.327 << "   "
		<< tbqsToEval(2)*197.327 << "   "
		<< tbqsToEval(3)*197.327 << "   "
		<< e*197.327 << "   " << eGiven*197.327 << "   "
		<< rhoB << "   " << rhoBGiven << "   "
		<< rhoS << "   " << rhoSGiven << "   "
		<< rhoQ << "   " << rhoQGiven << endl;
	*/
}
else
{
    e = (((rootfinder_parameters*)params)->eorEntSpline).eval(tbqsToEval);    //e, rhoB, rhoQ, rhoS contain the current point
    rhoB = (((rootfinder_parameters*)params)->rhoBSpline).eval(tbqsToEval);
    rhoQ = (((rootfinder_parameters*)params)->rhoQSpline).eval(tbqsToEval);
    rhoS = (((rootfinder_parameters*)params)->rhoSSpline).eval(tbqsToEval);
}



    gsl_vector_set(f, 0, (e - eGiven)); //f[0] contains (e(T,muB,muQ,muS) - eGiven)
    gsl_vector_set(f, 1, (rhoB - rhoBGiven)); //f[1] contains the (rhoB(T,muB,muQ,muS) - rhoBGiven)
    gsl_vector_set(f, 2, (rhoQ - rhoQGiven)); //f[2] contains the (rhoQ(T,muB,muQ,muS) - rhoQGiven)
    gsl_vector_set(f, 3, (rhoS - rhoSGiven)); //f[2] contains the (rho2(T,muB,muQ,muS) - rhoSGiven)

	/*std::cout << "Internal check(1): "
		<< gsl_vector_get(x,0) << "   " << gsl_vector_get(x,1) << "   "
		<< gsl_vector_get(x,3) << "   " << gsl_vector_get(x,2) << std::endl;
	std::cout << "Internal check(2): "
		<< e << "   " << rhoB << "   " << rhoS << "   " << rhoQ << std::endl;
	std::cout << "Internal check(3): "
		<< eGiven << "   " << rhoBGiven << "   " << rhoSGiven << "   " << rhoQGiven << std::endl;*/

    return GSL_SUCCESS;
}


bool eos::rootfinder4D(double e_or_s_Given, int e_or_s_mode, double rhoBGiven, double rhoSGiven, double rhoQGiven, double error, size_t steps)
{
	if ( VERBOSE > 5 ) std::cout << __PRETTY_FUNCTION__ << e_or_s_Given << "   " << e_or_s_mode << "   " << rhoBGiven << "   " << rhoSGiven << "   " << rhoQGiven << "   " << error << "   " << steps << std::endl;

	// Try the Delaunay interpolator first
	vector<double> result(4, 0.0);
	bool success = false;
	if ( e_or_s_mode==1 )
		success = e_delaunay.interpolate({e_or_s_Given*197.327, rhoBGiven, rhoSGiven, rhoQGiven}, result);
	else
		success = entr_delaunay.interpolate({e_or_s_Given, rhoBGiven, rhoSGiven, rhoQGiven}, result);

	{
		double phase_diagram_point[4] = {result[0], result[1], result[3], result[2]};	// NOTE: S <<-->> Q swapped!!!
		double densities_at_point[4];
		get_densities(phase_diagram_point, densities_at_point);
		cout << "Check solution:\n\t"
			<< phase_diagram_point[0] << "   "
			<< phase_diagram_point[1] << "   "
			<< phase_diagram_point[2] << "   "
			<< phase_diagram_point[3] << "\n\t"
			<< densities_at_point[0] << "   "
			<< densities_at_point[1] << "   "
			<< densities_at_point[2] << "   "
			<< densities_at_point[3] << "\n\t"
			<< e_or_s_Given*197.327 << "   "
			<< rhoBGiven << "   "
			<< rhoSGiven << "   "
			<< rhoQGiven << endl;
	}

	if ( success )
	{
		// set T, muB, muQ, muS
		tbqs( result[0]/197.327, result[1]/197.327, result[2]/197.327, result[3]/197.327 );
		return true;
	}

    //declare x = (T, muB, muQ, muS)
    gsl_vector *x = gsl_vector_alloc(4);

    gsl_vector_set(x, 0, T());
    gsl_vector_set(x, 1, muB());
    gsl_vector_set(x, 2, muQ());
    gsl_vector_set(x, 3, muS());
    /*gsl_vector_set(x, 0, 500.0/197.327);
    gsl_vector_set(x, 1, 0.0/197.327);
    gsl_vector_set(x, 2, 0.0/197.327);
    gsl_vector_set(x, 3, 0.0/197.327);*/

    //initialize the rootfinder equation to the correct variable quantities
    bool isEntropy = false;
    if(e_or_s_mode == 0) {
        isEntropy = true;
    }
    rootfinder_parameters p;
    if(isEntropy) {
        p.set(e_or_s_Given, rhoBGiven, rhoQGiven, rhoSGiven, entrSpline, bSpline, qSpline, sSpline);
    } else {
        p.set(e_or_s_Given, rhoBGiven, rhoQGiven, rhoSGiven, eSpline, bSpline, qSpline, sSpline);
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

	if ( e_or_s_mode == 1 and VERBOSE > 5 )
		std::cout << std::endl
			<< "=============================================="
			<< std::endl << "Input (e,B,Q,S): "
			<< e_or_s_Given*0.197327 << "   "
			<< rhoBGiven << "   "
			<< rhoQGiven << "   "
			<< rhoSGiven << std::endl;

    solver = gsl_multiroot_fsolver_alloc(TYPE, 4);
    gsl_multiroot_fsolver_set(solver, &f, x);

    int status;
    size_t iter = 0;

    do {
        ++iter;
        status = gsl_multiroot_fsolver_iterate(solver);


        if(VERBOSE > 5 && status) {

	if ( status == GSL_EBADFUNC && e_or_s_mode == 1 && VERBOSE > 5 )
		std::cout << "Error: something went to +/-Inf or NaN!" << std::endl;
	else if ( status == GSL_ENOPROG && e_or_s_mode == 1 && VERBOSE > 5 )
		std::cout << "Error: not making enough progress!" << std::endl;
            return 0;      //break if the rootfinder gets stuck
        }
        if(gsl_vector_get(solver->x, 0) < minT) {
	if ( e_or_s_mode == 1 && VERBOSE > 5 )
		std::cout << "Error: out-of-bounds (T < minT)!" << std::endl;
            return 0;
        } else if(gsl_vector_get(solver->x, 0) > maxT) {
	if ( e_or_s_mode == 1 && VERBOSE > 5 )
		std::cout << "Error: out-of-bounds (T > maxT)!" << std::endl;
            return 0;
        } else if (gsl_vector_get(solver->x, 1) < minMuB) {
	if ( e_or_s_mode == 1 && VERBOSE > 5 )
		std::cout << "Error: out-of-bounds (MuB < minMuB)!" << std::endl;
            return 0;
        } else if (gsl_vector_get(solver->x, 1) > maxMuB) {
	if ( e_or_s_mode == 1 && VERBOSE > 5 )
		std::cout << "Error: out-of-bounds (MuB > maxMuB)!" << std::endl;
            return 0;
        } else if (gsl_vector_get(solver->x, 2) < minMuQ) {     //break if the rootfinder goes out of bounds
	if ( e_or_s_mode == 1 && VERBOSE > 5 )
		std::cout << "Error: out-of-bounds (MuQ < minMuQ)!" << std::endl;
            return 0;
        } else if (gsl_vector_get(solver->x, 2) > maxMuQ) {
	if ( e_or_s_mode == 1 && VERBOSE > 5 )
		std::cout << "Error: out-of-bounds (MuQ > maxMuQ)!" << std::endl;
            return 0;
        } else if (gsl_vector_get(solver->x, 3) < minMuS) {
	if ( e_or_s_mode == 1 && VERBOSE > 5 )
		std::cout << "Error: out-of-bounds (MuS < minMuS)!" << std::endl;
            return 0;
        } else if (gsl_vector_get(solver->x, 3) > maxMuS) {
	if ( e_or_s_mode == 1 && VERBOSE > 5 )
		std::cout << "Error: out-of-bounds (MuS > maxMuS)!" << std::endl;
            return 0;
        }

        status = gsl_multiroot_test_residual(solver->f, error);

    } while (status == GSL_CONTINUE && iter < steps);


    bool found = true; //to return variable
    if(iter >= steps) {
        found = false;
    }


    if(found) {
        tbqs(gsl_vector_get(solver->x, 0), gsl_vector_get(solver->x, 1), gsl_vector_get(solver->x, 2), gsl_vector_get(solver->x, 3));    //set T, muB, muQ, muS
    }

	string output_status = ( found ) ? "FOUND" : "NOT FOUND";

	if ( e_or_s_mode == 1 && VERBOSE > 5 )
		std::cout << "Output (" << output_status << "): rootfinder4D at x = "
			<< 197.327*gsl_vector_get(solver->x, 0) << "   "
			<< 197.327*gsl_vector_get(solver->x, 1) << "   "
			<< 197.327*gsl_vector_get(solver->x, 2) << "   "
			<< 197.327*gsl_vector_get(solver->x, 3) << std::endl << std::endl;

    //memory deallocation
    gsl_multiroot_fsolver_free(solver);
    gsl_vector_free(x);
    return found;
}


