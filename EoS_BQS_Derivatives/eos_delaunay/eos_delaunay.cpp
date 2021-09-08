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
#include "delaunay.h"
#include "eos_delaunay.h"
#include "kdtree.h"
#include "point_in_simplex.h"

using namespace std;

eos_delaunay::eos_delaunay(string EoS_table_file)
{
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
	load_EoS_table(EoS_table_file, grid);

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

	cout << "Check minima and maxima:" << endl
		<< emin << "   " << emax << "   "
		<< bmin << "   " << bmax << "   "
		<< smin << "   " << smax << "   "
		<< qmin << "   " << qmax << endl;

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
	int emergency_count = 0;
	for (size_t iT   = 0; iT   < nT-1;   ++iT)
	for (size_t imub = 0; imub < nmub-1; ++imub)
	for (size_t imuq = 0; imuq < nmuq-1; ++imuq)
	for (size_t imus = 0; imus < nmus-1; ++imus)
	{
		std::array<double, 4> midpoint;
		for (size_t ii = 0; ii < 2; ++ii)
		for (size_t jj = 0; jj < 2; ++jj)
		for (size_t kk = 0; kk < 2; ++kk)
		for (size_t ll = 0; ll < 2; ++ll)
			std::transform( midpoint.begin(), midpoint.end(),
							grid[indexer( iT+ii, imub+jj, imuq+kk, imus+ll )].begin()+4,
							midpoint.begin(), std::plus<double>());

		std::transform( midpoint.begin(), midpoint.end(), midpoint.begin(),
							[](double & element){ return 0.0625*element; } );	//1/16

		midpoint_grid.push_back(midpoint);
		midpoint_inds.push_back( {iT, imub, imuq, imus} );
		if ( emergency_count++ <= 10 )
		{			
			cout << "Checking midpoint grid:" << endl;
			cout << "\t --> grid indices:";
			for (size_t ii = 0; ii < 2; ++ii)
			for (size_t jj = 0; jj < 2; ++jj)
			for (size_t kk = 0; kk < 2; ++kk)
			for (size_t ll = 0; ll < 2; ++ll)
				cout << " " << indexer( iT+ii, imub+jj, imuq+kk, imus+ll );
			cout << endl;
			cout << "\t --> phase diagram indices: "
					<< iT << "   " << imub << "   " << imuq << "   " << imus << endl;
			cout << "\t --> midpoint:";
			//for (auto elem : midpoint) cout << "   " << elem;
			cout << endl << endl;
		}

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
	for (size_t ii = 0; ii < grid.size(); ii++)
		std::copy_n( grid[ii].begin()+4, 4, density_points[ii].begin() );

	// set up kd-trees
	cout << "Setting up kd-trees...";
	static tree4d tree(std::begin(density_points), std::end(density_points));
	tree_ptr = &tree;

	
	cout << "Big check:" << endl;
	for ( const auto & midpoint : midpoint_grid )
	{
		for ( const auto & elem : midpoint )
			cout << "   " << elem;
		cout << endl;
	}
	
	static tree4d midpoint_tree(std::begin(midpoint_grid), std::end(midpoint_grid));
	midpoint_tree_ptr = &midpoint_tree;

	static tree4d unnormalized_midpoint_tree(std::begin(midpoint_unnormalized_grid),
											 std::end(midpoint_unnormalized_grid));
	unnormalized_midpoint_tree_ptr = &unnormalized_midpoint_tree;
	cout << "finished!\n";

	return;
}

void eos_delaunay::load_EoS_table(string path_to_file, vector<vector<double> > & grid)
{
	grid.clear();
	// read in file itself
	ifstream infile(path_to_file.c_str());
	if (infile.is_open())
	{
		size_t count = 0;
		string line;
		double dummy, Tin, muBin, muSin, muQin, bin, sin, qin, ein; 
		while ( getline (infile, line) )
		{
			count++;
			//if (count % 100 != 0) continue;

			istringstream iss(line);
			iss >> Tin >> muBin >> muQin >> muSin >> dummy >> dummy
				>> bin >> sin >> qin >> ein >> dummy;

			/*Tvec.push_back( Tin );
			muBvec.push_back( muBin );
			muSvec.push_back( muSin );
			muQvec.push_back( muQin );
			bvec.push_back( bin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) );		// 1/fm^3
			svec.push_back( sin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) );		// 1/fm^3
			qvec.push_back( qin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) );		// 1/fm^3
			evec.push_back( ein*Tin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) );	// MeV/fm^3
			*/

			grid.push_back( vector<double>({Tin, muBin, muQin, muSin,
											ein*Tin*Tin*Tin*Tin/(hbarc*hbarc*hbarc),
											bin*Tin*Tin*Tin/(hbarc*hbarc*hbarc),
											sin*Tin*Tin*Tin/(hbarc*hbarc*hbarc),
											qin*Tin*Tin*Tin/(hbarc*hbarc*hbarc)}) );

			if (count % 1000000 == 0)
			{
				cout << "Read in " << count << " lines." << endl;
				/*cout << "Most recently read in:" << endl
				<< Tin << "   " << muBin << "   " << muQin << "   " << muSin << "   " << 
				ein*Tin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) << "   " << 
				bin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) << "   " << 
				sin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) << "   " << 
				qin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) << endl;*/
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

void eos_delaunay::interpolate(const vector<double> & v0, vector<double> & result, bool verbose)
{
	if ( !interpolate_NMNmode(v0, result, verbose) )
	if ( !interpolate_NMNmode_v2(v0, result, false, verbose) )
	if ( !interpolate_NMNmode_v3(v0, result, verbose) )
		interpolate_NMNmode_v2(v0, result, true, verbose);
		//interpolate_NMNmode(v0, result);
}

// find containing simplex using nearest-midpoint-neighbor (NMN) method
bool eos_delaunay::interpolate_NMNmode(const vector<double> & v0, vector<double> & result, bool verbose)
{
//	cout << "Trying mode v1!" << endl;
	result.resize(4, 0.0);
	double e0 = v0[0], b0 = v0[1], s0 = v0[2], q0 = v0[3];

	// normalize first
	const double ne0 = (e0 - emin) / (emax - emin);
	const double nb0 = (b0 - bmin) / (bmax - bmin);
	const double ns0 = (s0 - smin) / (smax - smin);
	const double nq0 = (q0 - qmin) / (qmax - qmin);

	vector<double> nv0 = {ne0, nb0, ns0, nq0};

	// here is where we query the kd-tree for the nearest midpoint neighbor (NMN)
	size_t kdtree_nmn_index = 0;
	try
	{
		// point4d n not used; only need kdtree_nmn_index
		point4d n = midpoint_tree_ptr->nearest({ne0, nb0, ns0, nq0}, kdtree_nmn_index);
//		cout << "KD-Tree: NMN is " << n << endl;
//		cout << "KD-Tree: NMN distance: " << midpoint_tree_ptr->distance() << endl;
//		cout << "KD-Tree: NMN index is " << kdtree_nmn_index << endl;
//		cout << "KD-Tree: (T,muB,muQ,muS) indices of NMN are: "
//			<< midpoint_inds[kdtree_nmn_index][0] << ", "
//			<< midpoint_inds[kdtree_nmn_index][1] << ", "
//			<< midpoint_inds[kdtree_nmn_index][2] << ", "
//			<< midpoint_inds[kdtree_nmn_index][3] << endl;
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}

	// look up indices
	//const int iTNMN = Tinds[kdtree_nmn_index], imubNMN = mubinds[kdtree_nmn_index],
	//			imuqNMN = muqinds[kdtree_nmn_index], imusNMN = musinds[kdtree_nmn_index];
	const int iTNMN = midpoint_inds[kdtree_nmn_index][0];
	const int imubNMN = midpoint_inds[kdtree_nmn_index][1];
	const int imuqNMN = midpoint_inds[kdtree_nmn_index][2];
	const int imusNMN = midpoint_inds[kdtree_nmn_index][3];


	// select vertices in vicinity of NMN to triangulate
	int NMNvertex = 0;
	vector<vector<double> > vertices;

	// Qhull requires vertices as 1D vector
	vector<double> verticesFlat;

	// block for scope
	{
		int vertexcount = 0;
		for (int ii = 0; ii <= 1; ii++)
		for (int jj = 0; jj <= 1; jj++) // only need containing hypercube
		for (int kk = 0; kk <= 1; kk++) // vertices for the NMN method
		for (int ll = 0; ll <= 1; ll++)
		{
			// check that we're not going outside the grid
			if ( iTNMN+ii < nT && iTNMN+ii >= 0
				&& imubNMN+jj < nmub && imubNMN+jj >= 0
				&& imuqNMN+kk < nmuq && imuqNMN+kk >= 0
				&& imusNMN+ll < nmus && imusNMN+ll >= 0 )
			{
				if (ii==0 && jj==0 && kk==0 && ll==0)
					NMNvertex = vertexcount;	// identify NMN index below
				vertices.push_back( grid[indexer( iTNMN+ii, imubNMN+jj, imuqNMN+kk, imusNMN+ll )] );
				vertexcount++;
			}
		}

		size_t nVertices = vertices.size();
		/*if (nVertices < 6)	// this is how many Qhull needs
		{
				
			vertexcount = 0;
			vertices.clear();
			for (int ii = -1; ii <= 1; ii++)
			for (int jj = -1; jj <= 1; jj++) // only need containing hypercube
			for (int kk = -1; kk <= 1; kk++) // vertices for the NMN method
			for (int ll = -1; ll <= 1; ll++)
			{
				// check that we're not going outside the grid
				if ( iTNMN+ii < nT && iTNMN+ii >= 0
					&& imubNMN+jj < nmub && imubNMN+jj >= 0
					&& imuqNMN+kk < nmuq && imuqNMN+kk >= 0
					&& imusNMN+ll < nmus && imusNMN+ll >= 0 )
				{
					if (ii==0 && jj==0 && kk==0 && ll==0)
						NMNvertex = vertexcount;	// identify NMN index below
					vertices.push_back( grid[indexer( iTNMN+ii, imubNMN+jj, imuqNMN+kk, imusNMN+ll )] );
					vertexcount++;
				}
			}
			
			nVertices = vertices.size();
		}*/

		if (nVertices < 6) return false;	// just give up

		// flatten as efficiently as possible
		verticesFlat.resize(4*nVertices);	// dim == 4
		for (int ii = 0; ii < nVertices; ii++)
		{
			const vector<double> & vertex = vertices[ii];
			for (int jj = 0; jj < 4; jj++)
				verticesFlat[4*ii + jj] = vertex[jj+4];
		}
	
	}

	// Test the Delaunay part here
	// first get the triangulation
	vector<vector<size_t> > simplices;
	try
	{
		compute_delaunay(&verticesFlat[0], 4, verticesFlat.size() / 4, simplices);
	}
	catch (const std::exception& e)
	{
		std::cerr << '\n' << e.what() << '\n';
		std::cerr << __FUNCTION__ << ": Error occurred at "
				<< e0 << "   " << b0 << "   " << s0 << "   " << q0 << "\n";
		return false;
	}

	// =======================================================
	// triangulation is complete; now find containing simplex

	constexpr bool check_simplices = false;
	vector<bool> simplices_to_check(simplices.size(), false);

	// block to enforce local scope
	int iclosestsimplex = 0;
	{
		int isimplex = 0;
		double center_d2_min = 2.0;	// start with unrealistically large value (0 <= d2 <= 1)
		for ( const auto & simplex : simplices )
		{
			bool NMN_vertex_included_in_this_simplex = false;
			for ( const auto & vertex : simplex )
				if ( vertex == NMNvertex )
				{
					NMN_vertex_included_in_this_simplex = true;
					simplices_to_check[isimplex] = true;
					break;
				}
			
			// assume point must belong to simplex including NN, skip other simplices			
			if (check_simplices && !NMN_vertex_included_in_this_simplex)
			{
				isimplex++;
				continue;
			}

			// otherwise, compute simplex center and track squared distance to original point
			vector<double> center(4, 0.0);
			for ( const size_t vertex : simplex )
				std::transform( center.begin(), center.end(), vertices[vertex].begin()+4,
								center.begin(), std::plus<double>());


			// !!!!! N.B. - can remove this part and just multiply once !!!!!
			// !!!!! below by appropriate factors of 5					!!!!!
			// center is average of this simplex's vertices
			std::transform( center.begin(), center.end(), center.begin(),
							[](double & element){ return 0.2*element; } );
							// 0.2 == 1/(dim+1), dim == 4

			double d2loc = d2( center, nv0 );
			if ( d2loc < center_d2_min )
			{
				iclosestsimplex = isimplex;
				center_d2_min = d2loc;
			}
	
			isimplex++;
		}
	}

	// pass these to routine for locating point in simplex
	vector<vector<double> > simplexVertices(5);	// 5 == dim + 1, dim == 4

	// block for local scope
	{
		int ivertex = 0;
		for ( const auto & vertex : simplices[iclosestsimplex] )
			simplexVertices[ivertex++] = vector<double>( vertices[vertex].begin()+4,
									vertices[vertex].end() );
	}


	// try closest simplex first; otherwise loop through all simplices
	vector<double> point_lambda_in_simplex(5, 0.0);	// dim + 1 == 5
	bool foundPoint = point_is_in_simplex( simplexVertices, nv0, point_lambda_in_simplex, false );

	if (!foundPoint)        // loop over all simplices
	{
		int isimplex = 0;
		for ( auto & simplex : simplices )
		{
			if (check_simplices && !simplices_to_check[isimplex])
			{
				isimplex++;
				continue;       // skip simplices that don't need to be checked
			}

			// set simplex vertices
			simplexVertices.clear();
			for ( const auto & vertex : simplex )
				simplexVertices.push_back( vector<double>( vertices[vertex].begin()+4,
									vertices[vertex].end() ) );

			// check if point is in simplex; if so, return lambda coefficients and break
			if ( point_is_in_simplex( simplexVertices, nv0, point_lambda_in_simplex, false ) )
			{
				iclosestsimplex = isimplex;     // probably rename this
				foundPoint = true;
				break;
			}
			isimplex++;
		}
	}

	// finally, use the output lambda coefficients to get the interpolated values
	double T0 = 0.0, mub0 = 0.0, muq0 = 0.0, mus0 = 0.0;
	{
		int ivertex = 0;
		for ( const auto & vertex : simplices[iclosestsimplex] )
		{
			double lambda_coefficient = point_lambda_in_simplex[ivertex];
			T0   += lambda_coefficient * vertices[vertex][0];
			mub0 += lambda_coefficient * vertices[vertex][1];
			muq0 += lambda_coefficient * vertices[vertex][2];
			mus0 += lambda_coefficient * vertices[vertex][3];
			ivertex++;
		}
	}

	result[0] = T0;
	result[1] = mub0;
	result[2] = muq0;
	result[3] = mus0;

	return foundPoint;
}

// find containing simplex using nearest-midpoint-neighbor (NMN) method
bool eos_delaunay::interpolate_NMNmode_v2(
					const vector<double> & v0, vector<double> & result,
					bool expand_hypercube, bool verbose )
{
//	cout << "Trying mode v2!" << endl;
	result.resize(4, 0.0);
	double e0 = v0[0], b0 = v0[1], s0 = v0[2], q0 = v0[3];

	// normalize first
	const double ne0 = (e0 - emin) / (emax - emin);
	const double nb0 = (b0 - bmin) / (bmax - bmin);
	const double ns0 = (s0 - smin) / (smax - smin);
	const double nq0 = (q0 - qmin) / (qmax - qmin);

	vector<double> nv0 = {ne0, nb0, ns0, nq0};

	// here is where we query the kd-tree for the nearest midpoint neighbor (NMN)
	size_t kdtree_nmn_index = 0;
	try
	{
		// point4d n not used; only need kdtree_nmn_index
		point4d n = midpoint_tree_ptr->nearest({ne0, nb0, ns0, nq0}, kdtree_nmn_index);
		if (verbose)
		{
		size_t kdtree_nn_index = 0;
		point4d nn = tree_ptr->nearest({ne0, nb0, ns0, nq0}, kdtree_nn_index);
		cout << "KD-Tree: original point is "
				<< ne0 << "   " << nb0 << "   "
				<< ns0 << "   " << nq0 << endl;
		cout << "KD-Tree: NMN is " << n << endl;
		cout << "KD-Tree: NN is " << nn << endl;
		cout << "KD-Tree: NMN distance: " << midpoint_tree_ptr->distance() << endl;
		cout << "KD-Tree: NN distance: " << tree_ptr->distance() << endl;
		cout << "KD-Tree: NMN index is " << kdtree_nmn_index << endl;
		cout << "KD-Tree: NN index is " << kdtree_nn_index << endl;
		cout << "KD-Tree: (T,muB,muQ,muS) indices of NMN are: "
			<< midpoint_inds[kdtree_nmn_index][0] << ", "
			<< midpoint_inds[kdtree_nmn_index][1] << ", "
			<< midpoint_inds[kdtree_nmn_index][2] << ", "
			<< midpoint_inds[kdtree_nmn_index][3] << endl;
		cout << "KD-Tree: (T,muB,muQ,muS) indices of NN are: "
			<< Tinds[kdtree_nn_index] << ", "
			<< mubinds[kdtree_nn_index] << ", "
			<< muqinds[kdtree_nn_index] << ", "
			<< musinds[kdtree_nn_index] << endl;
		cout << "KD-Tree: vertices are:" << endl;
		for (int ii = 0; ii <= 1; ii++)
		for (int jj = 0; jj <= 1; jj++) // only need containing hypercube
		for (int kk = 0; kk <= 1; kk++) // vertices for the NMN method
		for (int ll = 0; ll <= 1; ll++)
		{
			if ( midpoint_inds[kdtree_nmn_index][0]+ii < nT
				&& midpoint_inds[kdtree_nmn_index][0]+ii >= 0
				&& midpoint_inds[kdtree_nmn_index][1]+jj < nmub
				&& midpoint_inds[kdtree_nmn_index][1]+jj >= 0
				&& midpoint_inds[kdtree_nmn_index][2]+kk < nmuq
				&& midpoint_inds[kdtree_nmn_index][2]+kk >= 0
				&& midpoint_inds[kdtree_nmn_index][3]+ll < nmus
				&& midpoint_inds[kdtree_nmn_index][3]+ll >= 0 )
			{
				const auto & cell = grid[indexer( midpoint_inds[kdtree_nmn_index][0]+ii,
								midpoint_inds[kdtree_nmn_index][1]+jj,
								midpoint_inds[kdtree_nmn_index][2]+kk,
								midpoint_inds[kdtree_nmn_index][3]+ll )];
				for (const auto & elem : cell)
					cout << "   " << elem;
				cout << endl;
			}
		}
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}

	// look up indices
	const int iTNMN = midpoint_inds[kdtree_nmn_index][0];
	const int imubNMN = midpoint_inds[kdtree_nmn_index][1];
	const int imuqNMN = midpoint_inds[kdtree_nmn_index][2];
	const int imusNMN = midpoint_inds[kdtree_nmn_index][3];


	// select vertices in vicinity of NMN to triangulate
	vector<vector<double> > vertices;

	int iclosestsimplex = 0;
	vector<vector<size_t> > simplices;
	vector<double> point_lambda_in_simplex(5, 0.0);	// dim + 1 == 5

	// big block around first attempt
	bool foundPoint = triangulate_and_locate_point( nv0,
						vector<int>({iTNMN, imubNMN, imuqNMN, imusNMN}),
						vertices, simplices, point_lambda_in_simplex, iclosestsimplex);


	// if we STILL have not found the containing simplex...
	if (expand_hypercube && !foundPoint)
	for (int iTshift = -1; iTshift <= -1; ++iTshift)
	for (int imubshift = -1; imubshift <= 1; ++imubshift)
	for (int imuqshift = -1; imuqshift <= 1; ++imuqshift)
	for (int imusshift = -1; imusshift <= 1; ++imusshift)
	{
		// the unshifted one was already tried above
		if (iTshift==0 && imubshift==0 && imuqshift==0 && imusshift==0) continue;

//		cout << "Trying " << iTNMN+iTshift << "   " << imubNMN+imubshift << "   "
//			<< imuqNMN+imuqshift << "   " << imusNMN+imusshift << ";   "
//			<< iTshift << "   " << imubshift << "   " << imuqshift << "   " << imusshift
//			<< endl;

		foundPoint = triangulate_and_locate_point( nv0,
						vector<int>({iTNMN+iTshift, imubNMN+imubshift,
									 imuqNMN+imuqshift, imusNMN+imusshift}),
						vertices, simplices, point_lambda_in_simplex, iclosestsimplex);

		if (foundPoint) goto end;
	}

	end:


	// finally, use the output lambda coefficients to get the interpolated values
	double T0 = 0.0, mub0 = 0.0, muq0 = 0.0, mus0 = 0.0;

	if (foundPoint)
	{
		int ivertex = 0;
		for ( const auto & vertex : simplices[iclosestsimplex] )
		{
			if ( vertex >= vertices.size() )
			{
				cerr << "Trying to access out of range!"<< endl;
				cerr << vertex << " vs. " << vertices.size() << endl;
			}
			double lambda_coefficient = point_lambda_in_simplex[ivertex];
			T0   += lambda_coefficient * vertices[vertex][0];
			mub0 += lambda_coefficient * vertices[vertex][1];
			muq0 += lambda_coefficient * vertices[vertex][2];
			mus0 += lambda_coefficient * vertices[vertex][3];
			ivertex++;
		}
	}
	

	result[0] = T0;
	result[1] = mub0;
	result[2] = muq0;
	result[3] = mus0;

	return foundPoint;
}



bool eos_delaunay::triangulate_and_locate_point(
					const vector<double> & nv0, const vector<int> & base,
					vector<vector<double> > & vertices, vector<vector<size_t> > & simplices,
					vector<double> & point_lambda_in_simplex, int & iclosestsimplex )
{
	const int iTbase = base[0], imubbase = base[1], imuqbase = base[2], imusbase = base[3];
	int NMNvertex = 0;

	vertices.clear();
	simplices.clear();

	// Qhull requires vertices as 1D vector
	vector<double> verticesFlat;

	// block for scope
	{
		int vertexcount = 0;
		for (int ii = 0; ii <= 1; ii++)
		for (int jj = 0; jj <= 1; jj++) // only need containing hypercube
		for (int kk = 0; kk <= 1; kk++) // vertices for the NMN method
		for (int ll = 0; ll <= 1; ll++)
//		for (int ii = -3; ii <= 4; ii++)
//		for (int jj = -3; jj <= 4; jj++) // only need containing hypercube
//		for (int kk = -3; kk <= 4; kk++) // vertices for the NMN method
//		for (int ll = -3; ll <= 4; ll++)
		{
			// check that we're not going outside the grid
			if ( iTbase+ii < nT && iTbase+ii >= 0
				&& imubbase+jj < nmub && imubbase+jj >= 0
				&& imuqbase+kk < nmuq && imuqbase+kk >= 0
				&& imusbase+ll < nmus && imusbase+ll >= 0 )
			{
				if (ii==0 && jj==0 && kk==0 && ll==0)
					NMNvertex = vertexcount;	// identify NMN index below
				vertices.push_back( grid[indexer( iTbase+ii, imubbase+jj, imuqbase+kk, imusbase+ll )] );
				vertexcount++;
			}
		}

		// flatten as efficiently as possible
		size_t nVertices = vertices.size();
		if (nVertices < 6) return false;

if (nVertices!=16)
{
	cout << "Not working with a true hypercube!  nVertices = " << nVertices << endl;
}
else
{
		refine_hypercube(vertices);	// just adds hypercube midpoint
		nVertices++;
}

		verticesFlat.resize(4*nVertices);	// dim == 4
		for (int ii = 0; ii < nVertices; ii++)
		{
			const vector<double> & vertex = vertices[ii];
			for (int jj = 0; jj < 4; jj++)
				verticesFlat[4*ii + jj] = vertex[jj+4];
		}
	
	}

	// Test the Delaunay part here
	// first get the triangulation
	//vector<vector<size_t> > simplices;
	try
	{
		compute_delaunay(&verticesFlat[0], 4, verticesFlat.size() / 4, simplices);
	}
	catch (const std::exception& e)
	{
		std::cerr << '\n' << e.what() << '\n';
		std::cerr << __FUNCTION__ << ": Error occurred at "
				<< nv0[0] << "   " << nv0[1] << "   " << nv0[2] << "   " << nv0[3] << "\n";
		return false;
	}

	// =======================================================
	// triangulation is complete; now find containing simplex

	constexpr bool check_simplices = false;
	vector<bool> simplices_to_check(simplices.size(), false);

	// block to enforce local scope
	iclosestsimplex = 0;
	{
		int isimplex = 0;
		double center_d2_min = 2.0;	// start with unrealistically large value (0 <= d2 <= 1)
		for ( const auto & simplex : simplices )
		{
			bool NMN_vertex_included_in_this_simplex = false;
			for ( const auto & vertex : simplex )
				if ( vertex == NMNvertex )
				{
					NMN_vertex_included_in_this_simplex = true;
					simplices_to_check[isimplex] = true;
					break;
				}
			
			// assume point must belong to simplex including NN, skip other simplices			
			if (check_simplices && !NMN_vertex_included_in_this_simplex)
			{
				isimplex++;
				continue;
			}

			// otherwise, compute simplex center and track squared distance to original point
			vector<double> center(4, 0.0);
			for ( const size_t vertex : simplex )
				std::transform( center.begin(), center.end(), vertices[vertex].begin()+4,
								center.begin(), std::plus<double>());


			// !!!!! N.B. - can remove this part and just multiply once !!!!!
			// !!!!! below by appropriate factors of 5					!!!!!
			// center is average of this simplex's vertices
			std::transform( center.begin(), center.end(), center.begin(),
							[](double & element){ return 0.2*element; } );
							// 0.2 == 1/(dim+1), dim == 4

			double d2loc = d2( center, nv0 );
			if ( d2loc < center_d2_min )
			{
				iclosestsimplex = isimplex;
				center_d2_min = d2loc;
			}
	
			isimplex++;
		}
	}

	// pass these to routine for locating point in simplex
	vector<vector<double> > simplexVertices(5);	// 5 == dim + 1, dim == 4

	// block for local scope
	{
		int ivertex = 0;
		for ( const auto & vertex : simplices[iclosestsimplex] )
			simplexVertices[ivertex++] = vector<double>( vertices[vertex].begin()+4,
									vertices[vertex].end() );
	}


	// try closest simplex first; otherwise loop through all simplices
	//vector<double> point_lambda_in_simplex(5, 0.0);	// dim + 1 == 5
	point_lambda_in_simplex.resize(5, 0.0);	// dim + 1 == 5
	bool foundPoint = point_is_in_simplex( simplexVertices, nv0, point_lambda_in_simplex, false );

//cout << "iclosestsimplex = " << iclosestsimplex << endl;

	if (!foundPoint)        // loop over all simplices
	{
		int isimplex = 0;
		for ( auto & simplex : simplices )
		{
//cout << "isimplex = " << isimplex << endl;

			if (check_simplices && !simplices_to_check[isimplex])
			{
				isimplex++;
				continue;       // skip simplices that don't need to be checked
			}

			// set simplex vertices
			simplexVertices.clear();
			for ( const auto & vertex : simplex )
				simplexVertices.push_back( vector<double>( vertices[vertex].begin()+4,
									vertices[vertex].end() ) );

			// check if point is in simplex; if so, return lambda coefficients and break
			if ( point_is_in_simplex( simplexVertices, nv0, point_lambda_in_simplex, false ) )
			{
				iclosestsimplex = isimplex;     // probably rename this
				foundPoint = true;
				break;
			}
			isimplex++;
		}
	}

	return foundPoint;
}

void eos_delaunay::refine_hypercube(vector<vector<double> > & hypercube)
{
	// first element assumed to represent lower corner of hypercube;
	// similarly for last element

	vector<double> lowerCorner(hypercube.front().begin(), hypercube.front().begin()+4);
	vector<double> upperCorner(hypercube.back().begin(), hypercube.back().begin()+4);
	vector<double> middle;
	std::transform( lowerCorner.begin(), lowerCorner.end(),
			upperCorner.begin(), std::back_inserter(middle),
			[](double x1, double x2) { return 0.5*(x1 + x2); });

	//for (int i0 = 0; i0 < 2; i0++)
	//for (int i1 = 0; i1 < 2; i1++)
	//for (int i2 = 0; i2 < 2; i2++)
	//for (int i3 = 0; i3 < 2; i3++)

	iter_swap(middle.begin() + 2, middle.begin() + 3);	// fix this eventually

//	cout << "lowerCorner:";
//	for (auto e : lowerCorner) cout << "   " << e;
//	cout << endl << "upperCorner:";
//	for (auto e : upperCorner) cout << "   " << e;
//	cout << endl << "middle:";
//	for (auto e : middle) cout << "   " << e;

	// just add midpoint for now	
	double densities_arr[4];
	get_densities(middle.data(), densities_arr);
	vector<double> densities(densities_arr, densities_arr+4);

	densities[0] = (densities[0] - emin)/(emax-emin);
	densities[1] = (densities[1] - bmin)/(bmax-bmin);
	densities[2] = (densities[2] - smin)/(smax-smin);
	densities[3] = (densities[3] - qmin)/(qmax-qmin);

	middle.insert( middle.end(), densities.begin(), densities.end() );

//	cout << endl << "middle:";
//	for (auto e : middle) cout << "   " << e;
//	cout << endl;

	hypercube.push_back( middle );

	return;
}

// find containing simplex using nearest-midpoint-neighbor (NMN) method
bool eos_delaunay::interpolate_NMNmode_v3(const vector<double> & v0, vector<double> & result, bool verbose)
{
//	cout << "Trying mode v3!" << endl;
	result.resize(4, 0.0);
	double e0 = v0[0], b0 = v0[1], s0 = v0[2], q0 = v0[3];

	// do not(!) normalize first
	const double ne0 = e0;
	const double nb0 = b0;
	const double ns0 = s0;
	const double nq0 = q0;

	vector<double> nv0 = {ne0, nb0, ns0, nq0};

	// here is where we query the kd-tree for the nearest midpoint neighbor (NMN)
	size_t kdtree_nmn_index = 0;
	try
	{
		// point4d n not used; only need kdtree_nmn_index
		point4d n = unnormalized_midpoint_tree_ptr->nearest(
					{ne0, nb0, ns0, nq0}, kdtree_nmn_index, false);	// false == log-distance mode
		if (verbose)
		{
		cout << "KD-Tree: original point is "
				<< ne0 << "   " << nb0 << "   "
				<< ns0 << "   " << nq0 << endl;
		cout << "KD-Tree: NMN is " << n << endl;
		cout << "KD-Tree: NMN distance: " << unnormalized_midpoint_tree_ptr->distance() << endl;
		cout << "KD-Tree: NMN index is " << kdtree_nmn_index << endl;
		cout << "KD-Tree: (T,muB,muQ,muS) indices of NMN are: "
			<< midpoint_inds[kdtree_nmn_index][0] << ", "
			<< midpoint_inds[kdtree_nmn_index][1] << ", "
			<< midpoint_inds[kdtree_nmn_index][2] << ", "
			<< midpoint_inds[kdtree_nmn_index][3] << endl;
		cout << "KD-Tree: vertices are:" << endl;
		for (int ii = 0; ii <= 1; ii++)
		for (int jj = 0; jj <= 1; jj++) // only need containing hypercube
		for (int kk = 0; kk <= 1; kk++) // vertices for the NMN method
		for (int ll = 0; ll <= 1; ll++)
		{
			if ( midpoint_inds[kdtree_nmn_index][0]+ii < nT
				&& midpoint_inds[kdtree_nmn_index][0]+ii >= 0
				&& midpoint_inds[kdtree_nmn_index][1]+jj < nmub
				&& midpoint_inds[kdtree_nmn_index][1]+jj >= 0
				&& midpoint_inds[kdtree_nmn_index][2]+kk < nmuq
				&& midpoint_inds[kdtree_nmn_index][2]+kk >= 0
				&& midpoint_inds[kdtree_nmn_index][3]+ll < nmus
				&& midpoint_inds[kdtree_nmn_index][3]+ll >= 0 )
			{
				const auto & cell
					= unnormalized_grid[indexer( midpoint_inds[kdtree_nmn_index][0]+ii,
								midpoint_inds[kdtree_nmn_index][1]+jj,
								midpoint_inds[kdtree_nmn_index][2]+kk,
								midpoint_inds[kdtree_nmn_index][3]+ll )];
				for (const auto & elem : cell)
					cout << "   " << elem;
				cout << endl;
			}
		}
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}

	// look up indices
	const int iTNMN = midpoint_inds[kdtree_nmn_index][0];
	const int imubNMN = midpoint_inds[kdtree_nmn_index][1];
	const int imuqNMN = midpoint_inds[kdtree_nmn_index][2];
	const int imusNMN = midpoint_inds[kdtree_nmn_index][3];


	// select vertices in vicinity of NMN to triangulate
	int NMNvertex = 0;
	vector<vector<double> > vertices;

	// Qhull requires vertices as 1D vector
	vector<double> verticesFlat;

	// block for scope
	{
		int vertexcount = 0;
		for (int ii = 0; ii <= 1; ii++)
		for (int jj = 0; jj <= 1; jj++) // only need containing hypercube
		for (int kk = 0; kk <= 1; kk++) // vertices for the NMN method
		for (int ll = 0; ll <= 1; ll++)
		{
			// check that we're not going outside the grid
			if ( iTNMN+ii < nT && iTNMN+ii >= 0
				&& imubNMN+jj < nmub && imubNMN+jj >= 0
				&& imuqNMN+kk < nmuq && imuqNMN+kk >= 0
				&& imusNMN+ll < nmus && imusNMN+ll >= 0 )
			{
				if (ii==0 && jj==0 && kk==0 && ll==0)
					NMNvertex = vertexcount;	// identify NMN index below
				size_t local_index = indexer( iTNMN+ii, imubNMN+jj, imuqNMN+kk, imusNMN+ll );
				vertices.push_back( unnormalized_grid[local_index] );
				vertexcount++;
			}
		}

		size_t nVertices = vertices.size();
		if (nVertices < 6) return false;	// just give up
/*if (nVertices!=16)
{
	cout << "Not working with a true hypercube!  nVertices = " << nVertices << endl;
}
else
{
		refine_hypercube(vertices);	// just adds hypercube midpoint
		nVertices++;
}*/

		// flatten as efficiently as possible
		verticesFlat.resize(4*nVertices);	// dim == 4
		for (int ii = 0; ii < nVertices; ii++)
		{
			const vector<double> & vertex = vertices[ii];
			for (int jj = 0; jj < 4; jj++)
				verticesFlat[4*ii + jj] = vertex[jj+4];
		}
	
	}

	// Test the Delaunay part here
	// first get the triangulation
	vector<vector<size_t> > simplices;
	try
	{
		compute_delaunay(&verticesFlat[0], 4, verticesFlat.size() / 4, simplices);
	}
	catch (const std::exception& e)
	{
		std::cerr << '\n' << e.what() << '\n';
		std::cerr << __FUNCTION__ << ": Error occurred at "
				<< e0 << "   " << b0 << "   " << s0 << "   " << q0 << "\n";
		return false;
	}

	// =======================================================
	// triangulation is complete; now find containing simplex

	constexpr bool check_simplices = false;
	vector<bool> simplices_to_check(simplices.size(), false);

	// block to enforce local scope
	int iclosestsimplex = 0;
	{
		int isimplex = 0;
		double center_d2_min = 2.0;	// start with unrealistically large value (0 <= d2 <= 1)
		for ( const auto & simplex : simplices )
		{
			bool NMN_vertex_included_in_this_simplex = false;
			for ( const auto & vertex : simplex )
				if ( vertex == NMNvertex )
				{
					NMN_vertex_included_in_this_simplex = true;
					simplices_to_check[isimplex] = true;
					break;
				}
			
			// assume point must belong to simplex including NN, skip other simplices			
			if (check_simplices && !NMN_vertex_included_in_this_simplex)
			{
				isimplex++;
				continue;
			}

			// otherwise, compute simplex center and track squared distance to original point
			vector<double> center(4, 0.0);
			for ( const size_t vertex : simplex )
				std::transform( center.begin(), center.end(), vertices[vertex].begin()+4,
								center.begin(), std::plus<double>());


			// !!!!! N.B. - can remove this part and just multiply once !!!!!
			// !!!!! below by appropriate factors of 5					!!!!!
			// center is average of this simplex's vertices
			std::transform( center.begin(), center.end(), center.begin(),
							[](double & element){ return 0.2*element; } );
							// 0.2 == 1/(dim+1), dim == 4

			double d2loc = d2( center, nv0 );
			if ( d2loc < center_d2_min )
			{
				iclosestsimplex = isimplex;
				center_d2_min = d2loc;
			}
	
			isimplex++;
		}
	}

	// pass these to routine for locating point in simplex
	vector<vector<double> > simplexVertices(5);	// 5 == dim + 1, dim == 4

	// block for local scope
	{
		int ivertex = 0;
		for ( const auto & vertex : simplices[iclosestsimplex] )
			simplexVertices[ivertex++] = vector<double>( vertices[vertex].begin()+4,
									vertices[vertex].end() );
	}


	// try closest simplex first; otherwise loop through all simplices
	vector<double> point_lambda_in_simplex(5, 0.0);	// dim + 1 == 5
	bool foundPoint = point_is_in_simplex( simplexVertices, nv0, point_lambda_in_simplex, false );

	if (!foundPoint)        // loop over all simplices
	{
		int isimplex = 0;
		for ( auto & simplex : simplices )
		{
			if (check_simplices && !simplices_to_check[isimplex])
			{
				isimplex++;
				continue;       // skip simplices that don't need to be checked
			}

			// set simplex vertices
			simplexVertices.clear();
			for ( const auto & vertex : simplex )
				simplexVertices.push_back( vector<double>( vertices[vertex].begin()+4,
									vertices[vertex].end() ) );

			// check if point is in simplex; if so, return lambda coefficients and break
			if ( point_is_in_simplex( simplexVertices, nv0, point_lambda_in_simplex, false ) )
			{
				iclosestsimplex = isimplex;     // probably rename this
				foundPoint = true;
				break;
			}
			isimplex++;
		}
	}

	// finally, use the output lambda coefficients to get the interpolated values
	double T0 = 0.0, mub0 = 0.0, muq0 = 0.0, mus0 = 0.0;
	{
		int ivertex = 0;
		for ( const auto & vertex : simplices[iclosestsimplex] )
		{
			double lambda_coefficient = point_lambda_in_simplex[ivertex];
			T0   += lambda_coefficient * vertices[vertex][0];
			mub0 += lambda_coefficient * vertices[vertex][1];
			muq0 += lambda_coefficient * vertices[vertex][2];
			mus0 += lambda_coefficient * vertices[vertex][3];
			ivertex++;
		}
	}

	result[0] = T0;
	result[1] = mub0;
	result[2] = muq0;
	result[3] = mus0;

	return foundPoint;
}



// find containing simplex using nearest-neighbor (NN) method
bool eos_delaunay::interpolate_NNmode(const vector<double> & v0, vector<double> & result, bool verbose)
{
	double e0 = v0[0], b0 = v0[1], s0 = v0[2], q0 = v0[3];

	// normalize first
	const double ne0 = (e0 - emin) / (emax - emin);
	const double nb0 = (b0 - bmin) / (bmax - bmin);
	const double ns0 = (s0 - smin) / (smax - smin);
	const double nq0 = (q0 - qmin) / (qmax - qmin);

	vector<double> nv0 = {ne0, nb0, ns0, nq0};

	// here is where we query the kd-tree for the nearest neighbor (NN)
	size_t kdtree_nn_index = 0;
	try
	{
		point4d n = tree_ptr->nearest({ne0, nb0, ns0, nq0}, kdtree_nn_index);
		//sw.Stop();
		//cout << "KD-Tree: Found nearest neighbor in " << setprecision(18)
		//		<< sw.printTime() << " s." << endl;
		//cout << "KD-Tree: Nearest neighbor is " << n << endl;
		//cout << "KD-Tree: Nearest neighbor index is " << kdtree_nn_index << endl;
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}

	// look up indices
	const int iTNN = Tinds[kdtree_nn_index], imubNN = mubinds[kdtree_nn_index],
				imuqNN = muqinds[kdtree_nn_index], imusNN = musinds[kdtree_nn_index];


	// select vertices in vicinity of NN to triangulate
	int NNvertex = 0;
	vector<vector<double> > vertices;

	// Qhull requires vertices as 1D vector
	vector<double> verticesFlat;

	// add block to constrain scope of unnecessary variables
	{
		int vertexcount = 0;
		for (int ii = -1; ii <= 1; ii++)
		for (int jj = -1; jj <= 1; jj++)
		for (int kk = -1; kk <= 1; kk++)
		for (int ll = -1; ll <= 1; ll++)
		{
			// check that we're not going outside the grid
			if ( iTNN+ii < nT && iTNN+ii >= 0
				&& imubNN+jj < nmub && imubNN+jj >= 0
				&& imuqNN+kk < nmuq && imuqNN+kk >= 0
				&& imusNN+ll < nmus && imusNN+ll >= 0 )
			{
				if (ii==0 && jj==0 && kk==0 && ll==0)
					NNvertex = vertexcount;	// identify NN index below
				vertices.push_back( grid[indexer( iTNN+ii, imubNN+jj, imuqNN+kk, imusNN+ll )] );
				vertexcount++;
			}
		}

		// flatten as efficiently as possible
		size_t nVertices = vertices.size();
		verticesFlat.resize(4*nVertices);	// dim == 4
		for (int ii = 0; ii < nVertices; ii++)
		{
			const vector<double> & vertex = vertices[ii];
			for (int jj = 0; jj < 4; jj++)
				verticesFlat[4*ii + jj] = vertex[jj+4];
		}
	}


	// Test the Delaunay part here
	// first get the triangulation
	vector<vector<size_t> > simplices;
	try
	{
		compute_delaunay(&verticesFlat[0], 4, verticesFlat.size() / 4, simplices);
		//cout << "Finished the Delaunay triangulation in " << sw.printTime() << " s." << endl;
	}
	catch (const std::exception& e)
	{
		std::cerr << '\n' << e.what() << '\n';
		std::cerr << __FUNCTION__ << ": Error occurred at "
				<< e0 << "   " << b0 << "   " << s0 << "   " << q0 << "\n";
		return false;
	}


	// =======================================================
	// triangulation is complete; now find containing simplex

	vector<bool> simplices_to_check(simplices.size(), false);

	// block to enforce local scope
	int iclosestsimplex = 0;
	{
		int isimplex = 0;
		double center_d2_min = 2.0;	// start with unrealistically large value (0 <= d2 <= 1)
		for ( const auto & simplex : simplices )
		{
			bool NN_vertex_included_in_this_simplex = false;
			for ( const auto & vertex : simplex )
				if ( vertex == NNvertex )
				{
					NN_vertex_included_in_this_simplex = true;
					simplices_to_check[isimplex] = true;
					break;
				}
			
			// assume point must belong to simplex including NN, skip other simplices			
			if (!NN_vertex_included_in_this_simplex)
			{
				isimplex++;
				continue;
			}

			// otherwise, compute simplex center and track squared distance to original point
			vector<double> center(4, 0.0);
			for ( const size_t vertex : simplex )
				std::transform( center.begin(), center.end(), vertices[vertex].begin()+4,
								center.begin(), std::plus<double>());


			// !!!!! N.B. - can remove this part and just multiple once !!!!!
			// !!!!! below by appropriate factors of 5					!!!!!
			// center is average of this simplex's vertices
			std::transform( center.begin(), center.end(), center.begin(),
							[](double & element){ return 0.2*element; } );
							// 0.2 == 1/(dim+1), dim == 4

			double d2loc = d2( center, nv0 );
			if ( d2loc < center_d2_min )
			{
				iclosestsimplex = isimplex;
				center_d2_min = d2loc;
			}
	
			isimplex++;
		}
	}

	// pass these to routine for locating point in simplex
	vector<vector<double> > simplexVertices(5);	// 5 == dim + 1, dim == 4

	// block for local scope
	{
		int ivertex = 0;
		for ( const auto & vertex : simplices[iclosestsimplex] )
			simplexVertices[ivertex++] = vector<double>( vertices[vertex].begin()+4,
														vertices[vertex].end() );
	}

	// locate the point in the simplex (assuming we know it's there;
	// currently no plan B for if point winds up outside this simplex)
	vector<double> point_lambda_in_simplex(5, 0.0);	// dim + 1 == 5
	bool foundPoint = point_is_in_simplex( simplexVertices, nv0, point_lambda_in_simplex, false );

	if (!foundPoint)        // loop over all simplices
	{
		//cout << "Did not find point in first simplex! Looping through all simplices:" << endl; 
		int isimplex = 0;
		for ( auto & simplex : simplices )
		{
			//cout << isimplex << ":" << endl;
			if (!simplices_to_check[isimplex])
			{
				isimplex++;
				continue;       // skip simplices that don't need to be checked
			}
			simplexVertices.clear();
			for ( const auto & vertex : simplex )
				simplexVertices.push_back( vector<double>( vertices[vertex].begin()+4,
									vertices[vertex].end() ) ); 
			if ( point_is_in_simplex( simplexVertices, nv0, point_lambda_in_simplex, false ) )
			{
				//cout << " found point in this simplex!" << endl;
				iclosestsimplex = isimplex;     // probably rename this
				foundPoint = true;
				break;
			}
			//else
			//	cout << " did not find point in this simplex!" << endl;
			isimplex++;
		}
	}

	// finally, use the output lambda coefficients to get the interpolated values
	double T0 = 0.0, mub0 = 0.0, muq0 = 0.0, mus0 = 0.0;
	{
		int ivertex = 0;
		for ( const auto & vertex : simplices[iclosestsimplex] )
		{
			double lambda_coefficient = point_lambda_in_simplex[ivertex];
			/*cout << "Check interpolation: " << ivertex << "   " << vertex << "   "
					<< point_lambda_in_simplex[ivertex] << "   "
					<< vertices[vertex][0] << "   " << vertices[vertex][1] << "   " 
					<< vertices[vertex][2] << "   " << vertices[vertex][3] << endl;*/
			T0   += lambda_coefficient * vertices[vertex][0];
			mub0 += lambda_coefficient * vertices[vertex][1];
			muq0 += lambda_coefficient * vertices[vertex][2];
			mus0 += lambda_coefficient * vertices[vertex][3];
			ivertex++;
		}
	}

	result.resize(4, 0.0);
	result[0] = T0;
	result[1] = mub0;
	result[2] = mus0;
	result[3] = muq0;

	return foundPoint;
}


