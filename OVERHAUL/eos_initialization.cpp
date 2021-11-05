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
#include "eos_delaunay/eos_delaunay.h"

#include "constants.h"

#include "rootfinder.h"

using namespace constants;

using std::vector;
using std::string;

constexpr bool use_exact = true;
constexpr bool accept_nearest_neighbor = false;
constexpr bool discard_unsolvable_charge_densities = false;

constexpr size_t STEPS = 1000000;
constexpr int VERBOSE = 0;
constexpr double TOLERANCE = 1e-12;


void EquationOfState::init()
{
  cout << "Attempting read in of EoS from "
        << quantity_file << " and " << deriv_file << endl;
  init( quantity_file, deriv_file );
}

void EquationOfState::init(string quantityFile, string derivFile)
{
	tbqsPosition.resize(4);

  // initialize things needed to use static C library
	cout << "Initializing EoS C library" << endl;
	initialize("/projects/jnorhos/BSQ/EoS_BQS_Derivatives/Coefficients_Parameters.dat");
  std::function<void(double[], double[])> f_eBSQ = get_eBSQ_densities;
  set_eBSQ_functional( f_eBSQ );
  std::function<void(double[], double[])> f_sBSQ = get_sBSQ_densities;
  set_sBSQ_functional( f_sBSQ );

  // load EoS tables, assess grid range
	std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
	init_grid_ranges_only(quantityFile, derivFile);

  // initialize Rootfinder ranges
  rootfinder.set_grid_ranges( minT, maxT, minMuB, maxMuB,
                              minMuS, maxMuS, minMuQ, maxMuQ );

  // initialize corresponding interpolator for each table
	cout << "Initialize Delaunay interpolators" << endl;
	e_delaunay.init(    quantityFile, 0 );	// 0 - energy density
	entr_delaunay.init( quantityFile, 1 );	// 1 - entropy density

	return;
}

void EquationOfState::init_grid_ranges_only(string quantityFile, string derivFile)
{
	if ( VERBOSE > 10 ) std::cout << "Now in " << __PRETTY_FUNCTION__ << std::endl;
    std::ifstream dataFile;
    dataFile.open(quantityFile);

    double tit, muBit, muQit, muSit, pit, entrit, bit, sit, qit, eit, cs2it;

    int count = 0;
    while ( dataFile >> tit >> muBit >> muQit >> muSit
                     >> pit >> entrit >> bit >> sit >> qit >> eit >> cs2it )
    {

		// Christopher Plumberg:
		// put T and mu_i in units of 1/fm
		tit   /= hbarc_MeVfm;
		muBit /= hbarc_MeVfm;
		muSit /= hbarc_MeVfm;
		muQit /= hbarc_MeVfm;

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