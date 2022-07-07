#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <lib.h>
#include "../include/eos_delaunay.h"
#include "../include/kdtree.h"

using namespace std;

const double eos_delaunay::hbarc = 197.3;
/*const size_t eos_delaunay::nT = 155,
			eos_delaunay::nmub = 37,
			eos_delaunay::nmus = 37,
			eos_delaunay::nmuq = 37;*/
const size_t eos_delaunay::nT = 241,
			eos_delaunay::nmub = 19,
			eos_delaunay::nmus = 19,
			eos_delaunay::nmuq = 19;
/*const size_t eos_delaunay::nT = 13,
			eos_delaunay::nmub = 9,
			eos_delaunay::nmus = 9,
			eos_delaunay::nmuq = 9;*/


//default constructor. This function exists to satisfy the compiler
//This function should never be called unless init is called directly afterward
eos_delaunay::eos_delaunay(){}

eos_delaunay::eos_delaunay(string EoS_table_file, int e_or_s)
{
	if ( e_or_s == 0 )
		cout << "\t --> initializing energy density" << endl;
	else
		cout << "\t --> initializing entropy density" << endl;
	init(EoS_table_file, e_or_s);
}

void eos_delaunay::init(string EoS_table_file, int e_or_s)
{
	using_e_or_s_mode = e_or_s;

	Tinds.resize(nT*nmub*nmuq*nmus);
	mubinds.resize(nT*nmub*nmuq*nmus);
	muqinds.resize(nT*nmub*nmuq*nmus);
	musinds.resize(nT*nmub*nmuq*nmus);

	size_t idx = 0;
	for (int iT = 0; iT < nT; iT++)
	for (int imub = 0; imub < nmub; imub++)
	for (int imuq = 0; imuq < nmuq; imuq++)
	for (int imus = 0; imus < nmus; imus++)
	{
		Tinds[idx] = iT;
		mubinds[idx] = imub;
		muqinds[idx] = imuq;
		musinds[idx] = imus;
		idx++;
	}

	// load EoS table
	load_EoS_table(EoS_table_file, grid, e_or_s);

	// check grid size; probably make this part of EoS table file's header info eventually...
	if ( grid.size() != nT*nmub*nmuq*nmus )
	{
		cerr << "Your chosen gridsizes do not match the input files!  Aborting!" << endl;
		exit(-1);
	}

	vector<double> Tvec(grid.size()), muBvec(grid.size()), muSvec(grid.size()), muQvec(grid.size());
	vector<double> evec(grid.size()), bvec(grid.size()), svec(grid.size()), qvec(grid.size());
	for (size_t igridcell = 0; igridcell < grid.size(); igridcell++)
	{
		const vector<double> & gridcell = grid[igridcell];
		Tvec[igridcell] = gridcell[0];
		muBvec[igridcell] = gridcell[1];
		muQvec[igridcell] = gridcell[2];	// note order of muQ, muS
		muSvec[igridcell] = gridcell[3];
		evec[igridcell] = gridcell[4];
		bvec[igridcell] = gridcell[5];
		svec[igridcell] = gridcell[6];		// note order of s, q
		qvec[igridcell] = gridcell[7];
	}

	// get density ranges and normalize
	emin = 0.0, emax = 0.0, bmin = 0.0, bmax = 0.0,
	smin = 0.0, smax = 0.0, qmin = 0.0, qmax = 0.0;
	get_min_and_max(evec, emin, emax, false);
	get_min_and_max(bvec, bmin, bmax, false);
	get_min_and_max(svec, smin, smax, false);
	get_min_and_max(qvec, qmin, qmax, false);

	// have one grid that we don't normalize for now
	unnormalized_grid = grid;

	// normalize grid points
	for ( vector<double> & gridcell : grid )
	{
		gridcell[4] = (gridcell[4] - emin) / ( emax - emin );
		gridcell[5] = (gridcell[5] - bmin) / ( bmax - bmin );
		gridcell[6] = (gridcell[6] - smin) / ( smax - smin );
		gridcell[7] = (gridcell[7] - qmin) / ( qmax - qmin );
	}


	// use midpoints as alternate way of find best simplex
	// "midpoints" are the average densities in the cell with lower corner at (iT,imu...)
	std::vector<std::array<double, 4> > midpoint_grid;
	for (size_t iT = 0; iT < nT-1; ++iT)
	for (size_t imub = 0; imub < nmub-1; ++imub)
	for (size_t imuq = 0; imuq < nmuq-1; ++imuq)
	for (size_t imus = 0; imus < nmus-1; ++imus)
	{
		std::array<double, 4> midpoint;
		midpoint.fill(0.0);
		for (size_t ii = 0; ii < 2; ++ii)
		for (size_t jj = 0; jj < 2; ++jj)
		for (size_t kk = 0; kk < 2; ++kk)
		for (size_t ll = 0; ll < 2; ++ll)
		{
			std::transform( midpoint.begin(), midpoint.end(),
							grid[indexer( iT+ii, imub+jj, imuq+kk, imus+ll )].begin()+4,
							midpoint.begin(), std::plus<double>());
		}

		std::transform( midpoint.begin(), midpoint.end(), midpoint.begin(),
							[](double & element){ return 0.0625*element; } );	//1/16

		midpoint_grid.push_back( midpoint );
		midpoint_inds.push_back( {iT, imub, imuq, imus} );
	}

	// use this for log-distance-based NMN method
	std::vector<std::array<double, 4> > midpoint_unnormalized_grid(midpoint_grid.size());
	{
		size_t iMidpoint = 0;
		for ( const std::array<double, 4> & midpoint : midpoint_grid )
		{
			midpoint_unnormalized_grid[iMidpoint][0] = emin + (emax-emin)*midpoint[0];
			midpoint_unnormalized_grid[iMidpoint][1] = bmin + (bmax-bmin)*midpoint[1];
			midpoint_unnormalized_grid[iMidpoint][2] = smin + (smax-smin)*midpoint[2];
			midpoint_unnormalized_grid[iMidpoint][3] = qmin + (qmax-qmin)*midpoint[3];
			iMidpoint++;
		}
	}


	// copy normalized densities from grid to separate vector (for kd-tree)
	// (needs to be vector of arrays to set up kd-tree correctly)
	std::vector<std::array<double, 4> > density_points(grid.size());
	std::vector<std::array<double, 4> > unnormalized_density_points(grid.size());
	for (size_t ii = 0; ii < grid.size(); ii++)
	{
		std::copy_n( grid[ii].begin()+4, 4, density_points[ii].begin() );
		std::copy_n( unnormalized_grid[ii].begin()+4, 4,
					 unnormalized_density_points[ii].begin() );
	}

	// set up kd-trees
	cout << "Setting up kd-trees...";
	if ( e_or_s == 0 )
	{
		static tree4d tree(std::begin(density_points), std::end(density_points));
		e_tree_ptr = &tree;
		
		static tree4d midpoint_tree(std::begin(midpoint_grid), std::end(midpoint_grid));
		e_midpoint_tree_ptr = &midpoint_tree;

		static tree4d unnormalized_tree(std::begin(unnormalized_density_points),
										std::end(unnormalized_density_points));
		e_unnormalized_tree_ptr = &unnormalized_tree;

		static tree4d unnormalized_midpoint_tree(std::begin(midpoint_unnormalized_grid),
												 std::end(midpoint_unnormalized_grid));
		e_unnormalized_midpoint_tree_ptr = &unnormalized_midpoint_tree;
	}
	else
	{
		static tree4d tree(std::begin(density_points), std::end(density_points));
		entr_tree_ptr = &tree;
		
		static tree4d midpoint_tree(std::begin(midpoint_grid), std::end(midpoint_grid));
		entr_midpoint_tree_ptr = &midpoint_tree;

		static tree4d unnormalized_tree(std::begin(unnormalized_density_points),
										std::end(unnormalized_density_points));
		entr_unnormalized_tree_ptr = &unnormalized_tree;

		static tree4d unnormalized_midpoint_tree(std::begin(midpoint_unnormalized_grid),
												 std::end(midpoint_unnormalized_grid));
		entr_unnormalized_midpoint_tree_ptr = &unnormalized_midpoint_tree;
	}
	cout << "finished!\n";

	return;
}

void eos_delaunay::load_EoS_table(string path_to_file, vector<vector<double> > & grid_in, int e_or_s)
{
	grid_in.clear();
	// read in file itself
	ifstream infile(path_to_file.c_str());
	if (infile.is_open())
	{
		size_t count = 0;
		string line;
		double dummy, Tin, muBin, muSin, muQin, bin, sin, qin, ein, entrin;
		if ( e_or_s == 0 )
		{
			while ( getline (infile, line) )
			{
				istringstream iss(line);
				iss >> Tin >> muBin >> muQin >> muSin >> dummy >> dummy
					>> bin >> sin >> qin >> ein >> dummy;

				grid_in.push_back( vector<double>({Tin, muBin, muQin, muSin,
												ein*Tin*Tin*Tin*Tin/(hbarc*hbarc*hbarc),
												bin*Tin*Tin*Tin/(hbarc*hbarc*hbarc),
												sin*Tin*Tin*Tin/(hbarc*hbarc*hbarc),
												qin*Tin*Tin*Tin/(hbarc*hbarc*hbarc)}) );

				if (++count % 1000000 == 0)
					cout << "Read in " << count << " lines." << endl;
			}
		}
		else
		{
			while ( getline (infile, line) )
			{
				istringstream iss(line);
				iss >> Tin >> muBin >> muQin >> muSin >> dummy >> entrin
					>> bin >> sin >> qin >> dummy >> dummy;

				grid_in.push_back( vector<double>({Tin, muBin, muQin, muSin,
												entrin*Tin*Tin*Tin/(hbarc*hbarc*hbarc),
												bin*Tin*Tin*Tin/(hbarc*hbarc*hbarc),
												sin*Tin*Tin*Tin/(hbarc*hbarc*hbarc),
												qin*Tin*Tin*Tin/(hbarc*hbarc*hbarc)}) );

				if (++count % 1000000 == 0)
					cout << "Read in " << count << " lines." << endl;
			}
		}
	}

	infile.close();
	return;
}

void eos_delaunay::get_min_and_max(vector<double> & v, double & minval, double & maxval, bool normalize)
{
	minval = *min_element(v.begin(), v.end());
	maxval = *max_element(v.begin(), v.end());
	if (normalize)
		std::transform( v.begin(), v.end(), v.begin(),
						[minval,maxval](double & element)
						{ return (element-minval)/(maxval-minval); } );
	return;
}


void eos_delaunay::get_NMN_coordinates(
					const vector<double> & v0, vector<double> & result,
					bool use_normalized)
{
	size_t kdtree_nmn_index = 0, kdtree_nn_index = 0;
	double nmn_dist = 2.0, nn_dist = 2.0;	// choose impossibly large separation
	try
	{
		if (use_normalized)
		{
			if ( using_e_or_s_mode == 0 )
			{
				midpoint_tree_ptr = e_midpoint_tree_ptr;
				tree_ptr          = e_tree_ptr;
			}
			else
			{
				midpoint_tree_ptr = entr_midpoint_tree_ptr;
				tree_ptr          = entr_tree_ptr;
			}
			
			point4d nmn = midpoint_tree_ptr->nearest
						( { (v0[0] - emin) / (emax - emin),
							(v0[1] - bmin) / (bmax - bmin),
							(v0[2] - smin) / (smax - smin),
							(v0[3] - qmin) / (qmax - qmin) }, kdtree_nmn_index );
			nmn_dist    = midpoint_tree_ptr->distance();
			point4d nn  = tree_ptr->nearest
						( { (v0[0] - emin) / (emax - emin),
							(v0[1] - bmin) / (bmax - bmin),
							(v0[2] - smin) / (smax - smin),
							(v0[3] - qmin) / (qmax - qmin) }, kdtree_nn_index );
			nn_dist     = tree_ptr->distance();
		}
		else	// using unnormalized grids
		{
			if ( using_e_or_s_mode == 0 )
			{
				unnormalized_midpoint_tree_ptr = e_unnormalized_midpoint_tree_ptr;
				unnormalized_tree_ptr          = e_unnormalized_tree_ptr;
			}
			else
			{
				unnormalized_midpoint_tree_ptr = entr_unnormalized_midpoint_tree_ptr;
				unnormalized_tree_ptr          = entr_unnormalized_tree_ptr;
			}
			
			point4d nmn = unnormalized_midpoint_tree_ptr->nearest
						( { v0[0], v0[1], v0[2], v0[3] }, kdtree_nmn_index );
			nmn_dist    = unnormalized_midpoint_tree_ptr->distance();
			point4d nn  = unnormalized_tree_ptr->nearest
						( { v0[0], v0[1], v0[2], v0[3] }, kdtree_nn_index );
			nn_dist     = unnormalized_tree_ptr->distance();
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}

	// approximate (T,muB,muQ,muS) coordinates of nearest midpoint neighbor
	if ( nmn_dist < nn_dist )
	{
		const vector<double> * vNMNptr = &grid[ indexer( midpoint_inds[kdtree_nmn_index] ) ];
		result.assign( vNMNptr->begin(), vNMNptr->begin()+4 );
	}
	else
	{
		const vector<double> * vNNptr = &grid[ kdtree_nn_index ];
		result.assign( vNNptr->begin(), vNNptr->end() );
	}
}


double eos_delaunay::normalized_d2(double a[], double b[])
{
	double na[4] = {(a[0]-emin)/(emax-emin), (a[1]-bmin)/(bmax-bmin),
					(a[2]-smin)/(smax-smin), (a[3]-qmin)/(qmax-qmin)};
	double nb[4] = {(b[0]-emin)/(emax-emin), (b[1]-bmin)/(bmax-bmin),
					(b[2]-smin)/(smax-smin), (b[3]-qmin)/(qmax-qmin)};
	return (  (na[0]-nb[0])*(na[0]-nb[0]) + (na[1]-nb[1])*(na[1]-nb[1])
			+ (na[2]-nb[2])*(na[2]-nb[2]) + (na[3]-nb[3])*(na[3]-nb[3]) );
}


double eos_delaunay::unnormalized_d2(double a[], double b[])
{
	return (  (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1])
			+ (a[2]-b[2])*(a[2]-b[2]) + (a[3]-b[3])*(a[3]-b[3]) );
}

double eos_delaunay::log_d2(double a[], double b[])
{
	double dist2 = 0;
	for (size_t iii = 0; iii < 4; ++iii)
	{
		double d = (a[iii]*b[iii]<=1e-20)? 1e10 : log(abs(a[iii]/b[iii]));
		dist2 += d * d;
	}
	return dist2;
}