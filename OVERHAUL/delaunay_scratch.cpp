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

using namespace constants;

using std::vector;
using std::string;

constexpr bool use_exact = true;
constexpr bool accept_nearest_neighbor = false;
constexpr bool discard_unsolvable_charge_densities = false;

constexpr size_t STEPS = 1000000;
constexpr int VERBOSE = 0;
constexpr double TOLERANCE = 1e-12;


void EquationOfState::delaunay_update_s
{

  // use NMN method to estimate where to start the rootfinder
	// ( returns estimates in units of MeV )
	vector<double> T_muB_muQ_muS_estimates;
	constexpr bool use_normalized_trees = true;
	if ( e_or_s_mode==1 )
		e_delaunay.get_NMN_coordinates(
					{e_or_s_Given*hbarc_MeVfm, rhoBGiven, rhoSGiven, rhoQGiven},
					T_muB_muQ_muS_estimates, use_normalized_trees );
	else
		entr_delaunay.get_NMN_coordinates(
					{e_or_s_Given, rhoBGiven, rhoSGiven, rhoQGiven},
					T_muB_muQ_muS_estimates, use_normalized_trees );




  if ( accept_nearest_neighbor )	// just return nearest neighbor instead
	{
		std::cout << "Entered to nearest neighbor tests" << std::endl;
		std::cout << "found = " << found << std::endl;

		// use unnormalized distances to estimate neighbor (reset from above)
		if ( e_or_s_mode==1 )
			e_delaunay.get_NMN_coordinates(
						{e_or_s_Given*hbarc_MeVfm, rhoBGiven, rhoSGiven, rhoQGiven},
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
			Tfinal   = gsl_vector_get(solver->x, 0)*hbarc_MeVfm;
			muBfinal = gsl_vector_get(solver->x, 1)*hbarc_MeVfm;
			muQfinal = gsl_vector_get(solver->x, 2)*hbarc_MeVfm;
			muSfinal = gsl_vector_get(solver->x, 3)*hbarc_MeVfm;
		}
		else
		{
			Tfinal   = previous_solver_step[0]*hbarc_MeVfm;
			muBfinal = previous_solver_step[1]*hbarc_MeVfm;
			muQfinal = previous_solver_step[2]*hbarc_MeVfm;
			muSfinal = previous_solver_step[3]*hbarc_MeVfm;
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
			tbqs( final_phase_diagram_point[0]/hbarc_MeVfm, final_phase_diagram_point[1]/hbarc_MeVfm,
				  final_phase_diagram_point[3]/hbarc_MeVfm, final_phase_diagram_point[2]/hbarc_MeVfm );
		cout << "final_phase_diagram_point: ";
		for (int iSol = 0; iSol <4; iSol++)
			cout << "   " << final_phase_diagram_point[iSol] / hbarc_MeVfm;
		cout << endl;

		}
		else if ( which_neighbor_closest == 1 )
		{
			tbqs( neighbor_estimate_point[0]/hbarc_MeVfm, neighbor_estimate_point[1]/hbarc_MeVfm,
				  neighbor_estimate_point[3]/hbarc_MeVfm, neighbor_estimate_point[2]/hbarc_MeVfm );
		cout << "neighbor_estimate_point!";
		for (int iSol = 0; iSol < 4; iSol++)
			cout << "   " << neighbor_estimate_point[iSol] / hbarc_MeVfm;
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

}