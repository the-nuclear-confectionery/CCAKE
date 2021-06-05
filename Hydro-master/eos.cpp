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

using std::vector;
using std::string;

using namespace SPLINTER;

// Compile:         gcc eos4D.cpp -c -I /usr/include/eigen3 -Lsplinter/build -lm -lgsl -lgslcblas -lstdc++ -lsplinter-3-0


//EoS constructor. Builds the splines of degree "degree" for each quantitiy and initializes the position at (30,0,0,0)
eos::eos(string quantityFile, string derivFile, int degree, bool using_HDF) : pSpline(4), entrSpline(4), bSpline(4), sSpline(4), qSpline(4), eSpline(4), cs2Spline(4), db2Spline(4), dq2Spline(4), ds2Spline(4), dt2Spline(4), dbdqSpline(4), dbdsSpline(4), dtdbSpline(4), dqdsSpline(4), dtdqSpline(4), dtdsSpline(4), tbqsPosition(4) {
    init(quantityFile, derivFile, degree, using_HDF);
}

//EoS default constructor. This function exists to satisfy the compiler
//This function should never be called unless init is called directly afterward
eos::eos() : pSpline(4), entrSpline(4), bSpline(4), sSpline(4), qSpline(4), eSpline(4), cs2Spline(4), db2Spline(4), dq2Spline(4), ds2Spline(4), dt2Spline(4), dbdqSpline(4), dbdsSpline(4), dtdbSpline(4), dqdsSpline(4), dtdqSpline(4), dtdsSpline(4), tbqsPosition(4) {}

void eos::init(string quantityFile, string derivFile, int degree, bool using_HDF)
{
	if ( using_HDF )
		init_with_hdf(quantityFile, derivFile, degree);
	else
		init_with_txt(quantityFile, derivFile, degree);

	return;
}


void eos::init_with_hdf(string quantityFile, string derivFile, int degree)
{
	const double hbarc = 197.327;

	Stopwatch sw;
	std::cout << "Loading HDF files...";
	sw.Start();
    vector<vector<double> > quantityData, derivData;
	read_in_hdf(quantityData, quantityFile);
	read_in_hdf(derivData, derivFile);
	sw.Stop();
	std::cout << "finished in " << sw.printTime() << " seconds!" << std::endl;

	std::cout << "Check dimensions: " << quantityData.size() << std::endl;
	std::cout << "Check dimensions: " << derivData.size() << std::endl;
	//if (true) exit(8);


	bool load_saved_files = true;

	if ( load_saved_files )
	{
		sw.Reset();
		sw.Start();
		std::cout << "Setting grid ranges...";
		bool initialize_ranges = true;
		for ( const auto & quantityRow : quantityData )
		{
			// put T and mu_i in units of 1/fm
			tit    = quantityRow[0]/hbarc;
			muBit  = quantityRow[1]/hbarc;
			muQit  = quantityRow[2]/hbarc;
			muSit  = quantityRow[3]/hbarc;		

		    if( initialize_ranges )
			{
		        minT   = tit;
		        maxT   = tit;
		        minMuB = muBit;
		        maxMuB = muBit;     //initialize eos range variables
		        minMuQ = muQit;
		        maxMuQ = muQit;
		        minMuS = muSit;
		        maxMuS = muSit;
				initialize_ranges = false;
		    }
			//if (count%100000==0) std::cout << "Read in line# " << count << std::endl;
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
		}
		sw.Stop();
		std::cout << "finished in " << sw.printTime() << " seconds!" << std::endl;


		sw.Reset();
		sw.Start();
		std::cout << "Loading all saved files:" << std::endl;
		/*std::cout << "Loading inputfiles/p.save..." << std::endl;
		psamples = DataTable("inputfiles/p.save");
		std::cout << "Loading inputfiles/entr.save..." << std::endl;
		entrsamples = DataTable("inputfiles/entr.save");
		std::cout << "Loading inputfiles/b.save..." << std::endl;
		bsamples = DataTable("inputfiles/b.save");
		std::cout << "Loading inputfiles/s.save..." << std::endl;
		ssamples = DataTable("inputfiles/s.save");
		std::cout << "Loading inputfiles/q.save..." << std::endl;
		qsamples = DataTable("inputfiles/q.save");
		std::cout << "Loading inputfiles/e.save..." << std::endl;
		esamples = DataTable("inputfiles/e.save");
		std::cout << "Loading inputfiles/cs2.save..." << std::endl;
		cs2samples = DataTable("inputfiles/cs2.save");
		
		std::cout << "Loading inputfiles/db2.save..." << std::endl;
		db2samples = DataTable("inputfiles/db2.save");
		std::cout << "Loading inputfiles/ds2.save..." << std::endl;
		ds2samples = DataTable("inputfiles/ds2.save");
		std::cout << "Loading inputfiles/dq2.save..." << std::endl;
		dq2samples = DataTable("inputfiles/dq2.save");
		std::cout << "Loading inputfiles/dt2.save..." << std::endl;
		dt2samples = DataTable("inputfiles/dt2.save");
		std::cout << "Loading inputfiles/dbds.save..." << std::endl;
		dbdssamples = DataTable("inputfiles/dbds.save");
		std::cout << "Loading inputfiles/dbdq.save..." << std::endl;
		dbdqsamples = DataTable("inputfiles/dbdq.save");
		std::cout << "Loading inputfiles/dqds.save..." << std::endl;
		dqdssamples = DataTable("inputfiles/dqds.save");
		std::cout << "Loading inputfiles/dtds.save..." << std::endl;
		dtdssamples = DataTable("inputfiles/dtds.save");
		std::cout << "Loading inputfiles/dtdq.save..." << std::endl;
		dtdqsamples = DataTable("inputfiles/dtdq.save");
		std::cout << "Loading inputfiles/dtdb.save..." << std::endl;
		dtdbsamples = DataTable("inputfiles/dtdb.save");*/


		std::cout << "Loading inputfiles/pSpline.save..." << std::endl;
		pSpline = BSpline("inputfiles/pSpline.save");
		std::cout << "Loading inputfiles/entrSpline.save..." << std::endl;
		entrSpline = BSpline("inputfiles/entrSpline.save");
		std::cout << "Loading inputfiles/bSpline.save..." << std::endl;
		bSpline = BSpline("inputfiles/bSpline.save");
		std::cout << "Loading inputfiles/sSpline.save..." << std::endl;
		sSpline = BSpline("inputfiles/sSpline.save");
		std::cout << "Loading inputfiles/qSpline.save..." << std::endl;
		qSpline = BSpline("inputfiles/qSpline.save");
		std::cout << "Loading inputfiles/eSpline.save..." << std::endl;
		eSpline = BSpline("inputfiles/eSpline.save");
		std::cout << "Loading inputfiles/cs2Spline.save..." << std::endl;
		cs2Spline = BSpline("inputfiles/cs2Spline.save");
		
		std::cout << "Loading inputfiles/db2Spline.save..." << std::endl;
		db2Spline = BSpline("inputfiles/db2Spline.save");
		std::cout << "Loading inputfiles/ds2Spline.save..." << std::endl;
		ds2Spline = BSpline("inputfiles/ds2Spline.save");
		std::cout << "Loading inputfiles/dq2Spline.save..." << std::endl;
		dq2Spline = BSpline("inputfiles/dq2Spline.save");
		std::cout << "Loading inputfiles/dt2Spline.save..." << std::endl;
		dt2Spline = BSpline("inputfiles/dt2Spline.save");
		std::cout << "Loading inputfiles/dbdsSpline.save..." << std::endl;
		dbdsSpline = BSpline("inputfiles/dbdsSpline.save");
		std::cout << "Loading inputfiles/dbdqSpline.save..." << std::endl;
		dbdqSpline = BSpline("inputfiles/dbdqSpline.save");
		std::cout << "Loading inputfiles/dqdsSpline.save..." << std::endl;
		dqdsSpline = BSpline("inputfiles/dqdsSpline.save");
		std::cout << "Loading inputfiles/dtdsSpline.save..." << std::endl;
		dtdsSpline = BSpline("inputfiles/dtdsSpline.save");
		std::cout << "Loading inputfiles/dtdqSpline.save..." << std::endl;
		dtdqSpline = BSpline("inputfiles/dtdqSpline.save");
		std::cout << "Loading inputfiles/dtdbSpline.save..." << std::endl;
		dtdbSpline = BSpline("inputfiles/dtdbSpline.save");

		sw.Stop();
		std::cout << "Finished loading all saved files in "
					<< sw.printTime() << " seconds!" << std::endl;
	}
	else
	{
	    DataTable psamples, entrsamples, bsamples, ssamples, qsamples, esamples, cs2samples;
	    DataTable db2samples, ds2samples, dq2samples, dt2samples, dbdssamples, dbdqsamples, dqdssamples, dtdssamples, dtdqsamples, dtdbsamples;

		double tit, muBit, muQit, muSit, pit, entrit, bit, sit, qit, eit, cs2it;
		double db2it, dq2it, ds2it, dt2it, dbdqit, dbdsit, dqdsit, dtdbit, dtdsit, dtdqit;
		vector<double> toAdd;

		Stopwatch sw_allocations, sw_addSample, sw_Total;

		long long count = 0;
		const long long nRows = quantityData.size();
		for ( long long iRow = 0; iRow < nRows; iRow++ )
		{
			sw_Total.Start();
			sw_allocations.Start();
			vector<double> & quantityRow = quantityData[iRow];
			vector<double> & derivRow = derivData[iRow];

			tit    = quantityRow[0];
			muBit  = quantityRow[1];
			muQit  = quantityRow[2];
			muSit  = quantityRow[3];
			pit    = quantityRow[4];
			entrit = quantityRow[5];
			bit    = quantityRow[6];
			sit    = quantityRow[7];
			qit    = quantityRow[8];
			eit    = quantityRow[9];
			cs2it  = quantityRow[10];
		
			db2it  = derivRow[4];
			dq2it  = derivRow[5];
			ds2it  = derivRow[6];
			dbdqit = derivRow[7];
			dbdsit = derivRow[8];
			dqdsit = derivRow[9];
			dtdbit = derivRow[10];
			dtdqit = derivRow[11];
			dtdsit = derivRow[12];
			dt2it  = derivRow[13];
			sw_allocations.Stop();
		

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

			// USE FM IN HYDRO
		    pit = pit*(tit*tit*tit*tit);
		    entrit = entrit*(tit*tit*tit);
		    bit = bit*(tit*tit*tit);
		    sit = sit*(tit*tit*tit);
		    qit = qit*(tit*tit*tit);
		    eit = eit*(tit*tit*tit*tit);
		

		    toAdd.push_back(tit);
		    toAdd.push_back(muBit);
		    toAdd.push_back(muQit);
		    toAdd.push_back(muSit);

			sw_addSample.Start();

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
			sw_addSample.Stop();

		    toAdd.clear();

			sw_Total.Stop();
			if (count%100000==0)
			{
				std::cout << "Spent " << sw_allocations.printTime()
							<< " seconds on allocations of total "
							<< sw_Total.printTime() << " seconds." << std::endl;
				std::cout << "Spent " << sw_addSample.printTime()
							<< " seconds on addSample of total "
							<< sw_Total.printTime() << " seconds." << std::endl;
			}
		}


		// try saving generated DataTables to files
		psamples.save("inputfiles/p.save");
		entrsamples.save("inputfiles/entr.save");
		bsamples.save("inputfiles/b.save");
		ssamples.save("inputfiles/s.save");
		qsamples.save("inputfiles/q.save");
		esamples.save("inputfiles/e.save");
		cs2samples.save("inputfiles/cs2.save");
		db2samples.save("inputfiles/db2.save");
		dq2samples.save("inputfiles/dq2.save");
		ds2samples.save("inputfiles/ds2.save");
		dbdqsamples.save("inputfiles/dbdq.save");
		dbdssamples.save("inputfiles/dbds.save");
		dqdssamples.save("inputfiles/dqds.save");
		dtdbsamples.save("inputfiles/dtdb.save");
		dtdqsamples.save("inputfiles/dtdq.save");
		dtdssamples.save("inputfiles/dtds.save");
		dt2samples.save("inputfiles/dt2.save");

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
	
	
		// save splines also
		pSpline.save("inputfiles/pSpline.save");
		entrSpline.save("inputfiles/entrSpline.save");
		bSpline.save("inputfiles/bSpline.save");
		sSpline.save("inputfiles/sSpline.save");
		qSpline.save("inputfiles/qSpline.save");
		eSpline.save("inputfiles/eSpline.save");
		cs2Spline.save("inputfiles/cs2Spline.save");
		db2Spline.save("inputfiles/db2Spline.save");
		dq2Spline.save("inputfiles/dq2Spline.save");
		ds2Spline.save("inputfiles/ds2Spline.save");
		dt2Spline.save("inputfiles/dt2Spline.save");
		dbdqSpline.save("inputfiles/dbdqSpline.save");
		dbdsSpline.save("inputfiles/dbdsSpline.save");
		dtdbSpline.save("inputfiles/dtdbSpline.save");
		dqdsSpline.save("inputfiles/dqdsSpline.save");
		dtdqSpline.save("inputfiles/dtdqSpline.save");
		dtdsSpline.save("inputfiles/dtdsSpline.save");
	}


	// initialize tbqsPosition to something...
	std::cout << "Initializing tbqsPosition...\n";
	for (int iTBQS = 0; iTBQS < 4; iTBQS++) tbqsPosition(iTBQS) = 1.0;

	std::cout << "Check TBQS: ";
	for (int iTBQS = 0; iTBQS < 4; iTBQS++) std::cout << tbqsPosition(iTBQS) << "   ";	
	std::cout << std::endl;

	std::cout << "All initializations finished!" << std::endl;


//if (true) exit(8);

    return;
}


void eos::init_with_txt(string quantityFile, string derivFile, int degree)
{
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
//	std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
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

/*void eos::set_tbqs(double setT, double setmuB, double setmuQ, double setmuS) {
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
}*/


double eos::T() {
	//std::cout << "Now in " << __PRETTY_FUNCTION__ << ": " << tbqsPosition(0) << std::endl;
    return tbqsPosition(0);
}

double eos::muB() {
	//std::cout << "Now in " << __PRETTY_FUNCTION__ << ": " << tbqsPosition(1) << std::endl;
    return tbqsPosition(1);
}

double eos::muQ() {
	//std::cout << "Now in " << __PRETTY_FUNCTION__ << ": " << tbqsPosition(2) << std::endl;
    return tbqsPosition(2);
}

double eos::muS() {
	//std::cout << "Now in " << __PRETTY_FUNCTION__ << ": " << tbqsPosition(3) << std::endl;
    return tbqsPosition(3);
}




double eos::p() {
    return pVal;
}

double eos::s() {
    return entrVal;
}

double eos::B() {
    return BVal;
}

double eos::S() {
    return SVal;
}

double eos::Q() {
    return QVal;
}

double eos::e() {
    return eVal;
}

double eos::cs2() {
    return cs2Val;
}


double eos::w() {
    return eVal + pVal;
}


/*double eos::p(double Tin, double muBin, double muQin, double muSin) {
	set_tbqs(setT, setmuB, setmuQ, setmuS);
    return pSpline.eval(tbqsPosition);
}

double eos::s(double Tin, double muBin, double muQin, double muSin) {
	set_tbqs(setT, setmuB, setmuQ, setmuS);
    return entrVal;
}

double eos::B(double Tin, double muBin, double muQin, double muSin) {
	set_tbqs(setT, setmuB, setmuQ, setmuS);
    return bSpline.eval(tbqsPosition);
}

double eos::S(double Tin, double muBin, double muQin, double muSin) {
	set_tbqs(setT, setmuB, setmuQ, setmuS);
    return sSpline.eval(tbqsPosition);
}

double eos::Q(double Tin, double muBin, double muQin, double muSin) {
    return qSpline.eval(tbqsPosition);
}

double eos::e(double Tin, double muBin, double muQin, double muSin) {
	set_tbqs(setT, setmuB, setmuQ, setmuS);
    return eSpline.eval(tbqsPosition);
}

double eos::cs2(double Tin, double muBin, double muQin, double muSin) {
	set_tbqs(setT, setmuB, setmuQ, setmuS);
    return cs2Spline.eval(tbqsPosition);
}


double eos::w(double Tin, double muBin, double muQin, double muSin) {
	set_tbqs(setT, setmuB, setmuQ, setmuS);
	pVal = pSpline.eval(tbqsPosition);
	eVal = eSpline.eval(tbqsPosition);
    return eVal + pVal;
}*/





double eos::dwds() {
    return T() + entrVal/dentr_dt() + BVal/dentr_dmub() + QVal/dentr_dmuq() + SVal/dentr_dmus();
}

double eos::dwdB() {
    return muB() + entrVal/db_dt() + BVal/db_dmub() + QVal/db_dmuq() + SVal/db_dmus();
}

double eos::dwdS() {
    return muS() + entrVal/ds_dt() + BVal/ds_dmub() + QVal/ds_dmuq() + SVal/ds_dmus();
}

double eos::dwdQ() {
    return muQ() + entrVal/dq_dt() + BVal/dq_dmub() + QVal/dq_dmuq() + SVal/dq_dmus();
}

double eos::dentr_dt() {
    return calc_term_1();
}

double eos::dentr_dmub() {
    return calc_term_2("b");
}

double eos::dentr_dmuq() {
    return calc_term_2("q");
}

double eos::dentr_dmus() {
    return calc_term_2("s");
}

double eos::db_dt() {
    return calc_term_3("b");
}

double eos::db_dmub() {
    return calc_term_4("b","b");
}

double eos::db_dmuq() {
    return calc_term_4("b","q");
}

double eos::db_dmus() {
    return calc_term_4("b","s");
}

double eos::ds_dt() {
    return calc_term_3("s");
}

double eos::ds_dmub() {
    return calc_term_4("s","b");
}

double eos::ds_dmuq() {
    return calc_term_4("s","q");
}

double eos::ds_dmus() {
    return calc_term_4("s","s");
}

double eos::dq_dt() {
    return calc_term_3("q");
}

double eos::dq_dmub() {
    return calc_term_4("q","b");
}

double eos::dq_dmuq() {
    return calc_term_4("q","q");
}

double eos::dq_dmus() {
    return calc_term_4("q","s");
}

double eos::calc_term_1() {
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
//	std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
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
//	std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
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
//	std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
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
//	std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    gsl_permutation *p = gsl_permutation_alloc(3);
    int s;

    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(m, p, &s);

    // Compute the  inverse of the LU decomposition
    gsl_matrix *minv = gsl_matrix_alloc(3, 3);
    gsl_linalg_LU_invert(m, p, minv);

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
//	std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    Aout=w()-entrVal*dwds();

    return Aout;
}

double eos::cs2out(double Tt) {  //return cs2 given t and mu's=0
//	std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    tbqs(Tt, 0.0, 0.0, 0.0);
    return cs2Val;
}

double eos::cs2out(double Tt, double muBin, double muQin, double muSin) {  //return cs2 given t and mu's
//	std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    tbqs(Tt, muBin, muQin, muSin);
    return cs2Val;
}

double eos::wfz(double Tt) {   // return e + p for tbqs
//	std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    tbqs(Tt, 0.0, 0.0, 0.0);
    return eVal + pVal;
}

double eos::wfz(double Tt, double muBin, double muQin, double muSin) {   // return e + p for tbqs
//	std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    tbqs(Tt, muBin, muQin, muSin);
    return eVal + pVal;
}

bool eos::update_s(double sin) { //update the t position (mu=0) based on input. Returns 1 if found, returns 0 if failed
//	std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return update_s(sin, 0.0, 0.0, 0.0);
}

bool eos::update_s(double sin, double Bin, double Sin, double Qin) { //update the t and mu position based on input. Returns 1 if found, returns 0 if failed
//	std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
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
//	std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    return s_out(ein, 0.0, 0.0, 0.0);
}

double eos::s_out(double ein, double Bin, double Sin, double Qin) {   //update the t and mu position based on input. Returns entropy if found, returns -1 if failed
//	std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
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
//	std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    tbqs(Tt, 0, 0, 0);
    return entrVal;
}


// UNCOMMENTED BY C. PLUMBERG
void eos::eosin(std::string type) {}
double eos::A() {
    return 0;
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
    e = (((rootfinder_parameters*)params)->eorEntSpline).eval(tbqsToEval);    //e, rhoB, rhoQ, rhoS contain the current point
    rhoB = (((rootfinder_parameters*)params)->rhoBSpline).eval(tbqsToEval);
    rhoQ = (((rootfinder_parameters*)params)->rhoQSpline).eval(tbqsToEval);
    rhoS = (((rootfinder_parameters*)params)->rhoSSpline).eval(tbqsToEval);




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



bool eos::rootfinder4D(double e_or_s_Given, int e_or_s_mode, double rhoBGiven, double rhoSGiven, double rhoQGiven, double error, size_t steps) {
//	std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;

    //declare x = (T, muB, muS)
    gsl_vector *x = gsl_vector_alloc(4);
    gsl_vector_set(x, 0, T());
    gsl_vector_set(x, 1, muB());
    gsl_vector_set(x, 2, muQ());
    gsl_vector_set(x, 3, muS());

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

	if ( e_or_s_mode == 1 )
		std::cout << std::endl
			<< "=============================================="
			<< std::endl << "Input (e,B,Q,S): "
			<< e_or_s_Given*0.19733 << "   "
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


        if(status) {

	if ( status == GSL_EBADFUNC && e_or_s_mode == 1 )
		std::cout << "Error: something went to +/-Inf or NaN!" << std::endl;
	else if ( status == GSL_ENOPROG && e_or_s_mode == 1 )
		std::cout << "Error: not making enough progress!" << std::endl;
            return 0;      //break if the rootfinder gets stuck
        }
        if(gsl_vector_get(solver->x, 0) < minT) {
	if ( e_or_s_mode == 1 )
		std::cout << "Error: out-of-bounds (T < minT)!" << std::endl;
            return 0;
        } else if(gsl_vector_get(solver->x, 0) > maxT) {
	if ( e_or_s_mode == 1 )
		std::cout << "Error: out-of-bounds (T > maxT)!" << std::endl;
            return 0;
        } else if (gsl_vector_get(solver->x, 1) < minMuB) {
	if ( e_or_s_mode == 1 )
		std::cout << "Error: out-of-bounds (MuB < minMuB)!" << std::endl;
            return 0;
        } else if (gsl_vector_get(solver->x, 1) > maxMuB) {
	if ( e_or_s_mode == 1 )
		std::cout << "Error: out-of-bounds (MuB > maxMuB)!" << std::endl;
            return 0;
        } else if (gsl_vector_get(solver->x, 2) < minMuQ) {     //break if the rootfinder goes out of bounds
	if ( e_or_s_mode == 1 )
		std::cout << "Error: out-of-bounds (MuQ < minMuQ)!" << std::endl;
            return 0;
        } else if (gsl_vector_get(solver->x, 2) > maxMuQ) {
	if ( e_or_s_mode == 1 )
		std::cout << "Error: out-of-bounds (MuQ > maxMuQ)!" << std::endl;
            return 0;
        } else if (gsl_vector_get(solver->x, 3) < minMuS) {
	if ( e_or_s_mode == 1 )
		std::cout << "Error: out-of-bounds (MuS < minMuS)!" << std::endl;
            return 0;
        } else if (gsl_vector_get(solver->x, 3) > maxMuS) {
	if ( e_or_s_mode == 1 )
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

	if ( e_or_s_mode == 1 )
		std::cout << "Output (" << output_status << "): rootfinder4D at x = "
			<< 197.33*gsl_vector_get(solver->x, 0) << "   "
			<< 197.33*gsl_vector_get(solver->x, 1) << "   "
			<< 197.33*gsl_vector_get(solver->x, 2) << "   "
			<< 197.33*gsl_vector_get(solver->x, 3) << std::endl << std::endl;

    //memory deallocation
    gsl_multiroot_fsolver_free(solver);
    gsl_vector_free(x);
    return found;
}







//struct to pass the target quantities into the rootfinder function
struct quant_rootfinder_parameters {
    double tGiven;
    double muBGiven;
    double muQGiven;
    double muSGiven;
    double quantGiven;          //the value we are looking for in the desired quantity
    int whichIndep;
    BSpline pSpline;        //the spine of the correct quantity
    BSpline entrSpline;
    BSpline eSpline;
    string quantity;
    quant_rootfinder_parameters();
    quant_rootfinder_parameters(string setQuantity, int setWhichIndep, double setQuantGiven, double setT, double setmuB, double setmuQ, double setmuS, BSpline setPSpline, BSpline setEntrSpline, BSpline setESpline);
public:
    void set(string setQuantity, int setWhichIndep, double setQuantGiven, double setT, double setmuB, double setmuQ, double setmuS, BSpline setPSpline, BSpline setEntrSpline, BSpline setESpline);
};
//Default constructor
quant_rootfinder_parameters::quant_rootfinder_parameters() : pSpline(4), entrSpline(4), eSpline(4) {}
//constructor which initializes all struct variables
quant_rootfinder_parameters::quant_rootfinder_parameters(string setQuantity, int setWhichIndep, double setQuantGiven, double setT, double setmuB, double setmuQ, double setmuS, BSpline setPSpline, BSpline setEntrSpline, BSpline setESpline) : pSpline(4), entrSpline(4), eSpline(4)
{
    whichIndep = setWhichIndep;
    quantity = setQuantity;
    quantGiven = setQuantGiven;
    pSpline = setPSpline;
    entrSpline = setEntrSpline;
    eSpline = setESpline;
    tGiven = setT;
    muBGiven = setmuB;
    muQGiven = setmuQ;
    muSGiven = setmuS;
}
void quant_rootfinder_parameters::set(string setQuantity, int setWhichIndep, double setQuantGiven, double setT, double setmuB, double setmuQ, double setmuS, BSpline setPSpline, BSpline setEntrSpline, BSpline setESpline) {
    whichIndep = setWhichIndep;
    quantity = setQuantity;
    quantGiven = setQuantGiven;
    pSpline = setPSpline;
    entrSpline = setEntrSpline;
    eSpline = setESpline;
    tGiven = setT;
    muBGiven = setmuB;
    muQGiven = setmuQ;
    muSGiven = setmuS;
}

//helper function for the rootfinder. It provides the correct difference of quantity from the target
int quant_rootfinder_f(const gsl_vector *x, void *params, gsl_vector *f);
int quant_rootfinder_f(const gsl_vector *x, void *params, gsl_vector *f) {
    int whichIndep = ((quant_rootfinder_parameters*)params)->whichIndep;
    //x contains the next (T, muB, muS) coordinate to test
    DenseVector tbqsToEval(4);
    if(whichIndep == 1) {
        tbqsToEval(0) = gsl_vector_get(x,0);
        tbqsToEval(1) = ((quant_rootfinder_parameters*)params)->muBGiven;      //convert x into densevector so it can be a BSpline evaluation point
        tbqsToEval(2) = ((quant_rootfinder_parameters*)params)->muQGiven;
        tbqsToEval(3) = ((quant_rootfinder_parameters*)params)->muSGiven;
    } else if(whichIndep == 2) {
        tbqsToEval(0) = ((quant_rootfinder_parameters*)params)->tGiven;
        tbqsToEval(1) = gsl_vector_get(x,0);      //convert x into densevector so it can be a BSpline evaluation point
        tbqsToEval(2) = ((quant_rootfinder_parameters*)params)->muQGiven;
        tbqsToEval(3) = ((quant_rootfinder_parameters*)params)->muSGiven;
    } else if(whichIndep == 3) {
        tbqsToEval(0) = ((quant_rootfinder_parameters*)params)->tGiven;
        tbqsToEval(1) = ((quant_rootfinder_parameters*)params)->muBGiven;      //convert x into densevector so it can be a BSpline evaluation point
        tbqsToEval(2) = gsl_vector_get(x,0);
        tbqsToEval(3) = ((quant_rootfinder_parameters*)params)->muSGiven;
    } else {
        tbqsToEval(0) = ((quant_rootfinder_parameters*)params)->tGiven;
        tbqsToEval(1) = ((quant_rootfinder_parameters*)params)->muBGiven;      //convert x into densevector so it can be a BSpline evaluation point
        tbqsToEval(2) = ((quant_rootfinder_parameters*)params)->muQGiven;
        tbqsToEval(3) = gsl_vector_get(x,0);
    }



    double quantGiven, quant;
    if(((quant_rootfinder_parameters*)params)->quantity == "e") {
        quantGiven = ((quant_rootfinder_parameters*)params)->quantGiven;
        quant = (((quant_rootfinder_parameters*)params)->eSpline).eval(tbqsToEval);
    } else if(((quant_rootfinder_parameters*)params)->quantity == "p") {
        quantGiven = ((quant_rootfinder_parameters*)params)->quantGiven;
        quant = (((quant_rootfinder_parameters*)params)->pSpline).eval(tbqsToEval);
    } else if(((quant_rootfinder_parameters*)params)->quantity == "entr") {
        quantGiven = ((quant_rootfinder_parameters*)params)->quantGiven;
        quant = (((quant_rootfinder_parameters*)params)->entrSpline).eval(tbqsToEval);
    } else if(((quant_rootfinder_parameters*)params)->quantity == "gibbs") {
        quantGiven = ((quant_rootfinder_parameters*)params)->quantGiven;
        quant = (((quant_rootfinder_parameters*)params)->eSpline).eval(tbqsToEval);
        quant += (((quant_rootfinder_parameters*)params)->pSpline).eval(tbqsToEval);
        quant -= ((((quant_rootfinder_parameters*)params)->entrSpline).eval(tbqsToEval))*tbqsToEval(0);
    }

    gsl_vector_set(f, 0, (quant - quantGiven));

    return GSL_SUCCESS;
}

bool eos::quant_rootfinder4D(double quantGiven, string quantType, int whichIndep, double error, size_t steps) {
//	std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;

    //declare x = T
    gsl_vector *x = gsl_vector_alloc(1);
    quant_rootfinder_parameters p;
    if(whichIndep == 1) {
        gsl_vector_set(x, 0, T());
    } else if(whichIndep == 2) {
        gsl_vector_set(x, 0, muB());
    } else if(whichIndep == 3) {

        gsl_vector_set(x, 0, muQ());
    } else if(whichIndep == 4) {
        gsl_vector_set(x, 0, muS());
    } else {
        std::cout << "Please select a gridded quantity to vary during search" << std::endl;
        std::cout << "t = 1,  muB = 2,  muQ = 3,  muS = 4" << std::endl;
        return 0;
    }
    p.set(quantType, whichIndep, quantGiven, T(), muB(), muQ(), muS(), pSpline, entrSpline, eSpline);

    //initialize multiroot solver
    const gsl_multiroot_fsolver_type *type;
    gsl_multiroot_fsolver *solver;
    gsl_multiroot_function f;
    f.f = &quant_rootfinder_f;
    f.n = 1;
    f.params = &p;

    //type options: dnewton, hybrids, hybrid, broyden
    //the dnewton and broyden methods crash the program with a bad guess. Hybrid and hybrids do not
    //dnewton works the fastest on guesses that it can successfully reach. Hybrids is faster than hybrid
    //I am choosing hybrids because it is the most accurate. Hybrids is slower than dnewton, but more reliable
    type = gsl_multiroot_fsolver_hybrids;
    solver = gsl_multiroot_fsolver_alloc(type, 1);
    gsl_multiroot_fsolver_set(solver, &f, x);

    int status;
    size_t iter = 0;

    do {
        ++iter;
        status = gsl_multiroot_fsolver_iterate(solver);

        if(status) {
            return 0;      //break if the rootfinder gets stuck
        }
        if(whichIndep == 1) {
            if(gsl_vector_get(solver->x, 0) < minT) {
                return 0;
            } else if(gsl_vector_get(solver->x, 0) > maxT) {
                return 0;
            }
        } else if(whichIndep == 2) {
            if(gsl_vector_get(solver->x, 0) < minMuB) {
                return 0;
            } else if(gsl_vector_get(solver->x, 0) > maxMuB) {
                return 0;
            }
        } else if(whichIndep == 3) {
            if(gsl_vector_get(solver->x, 0) < minMuQ) {
                return 0;
            } else if(gsl_vector_get(solver->x, 0) > maxMuQ) {
                return 0;
            }
        } else if(whichIndep == 4) {
            if(gsl_vector_get(solver->x, 0) < minMuS) {
                return 0;
            } else if(gsl_vector_get(solver->x, 0) > maxMuS) {
                return 0;
            }
        }


        status = gsl_multiroot_test_residual(solver->f, error);

    } while (status == GSL_CONTINUE && iter < steps);


    bool found = true; //to return variable
    if(iter >= steps) {
        found = false;
    }


    if(found) {
        tbqs(gsl_vector_get(solver->x, 0), muB(), muQ(), muS());    //set T, muB, muQ, muS
    }

    //memory deallocation
    gsl_multiroot_fsolver_free(solver);
    gsl_vector_free(x);
    return found;
}
