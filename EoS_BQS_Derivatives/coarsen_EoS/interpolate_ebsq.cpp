#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string>

#include "delaunay.h"
#include "kdtree.h"
#include "point_in_simplex.h"
#include "Stopwatch.h"

using namespace std;

constexpr int mode = 1;
constexpr double EPS = 1e-25;
const double hbarc = 197.327;

constexpr size_t nT = 155, nmub = 37, nmus = 37, nmuq = 37;

inline size_t indexer( const int iT, const int imub, const int imus, const int imuq )
{
	// mus varies faster than muq!!!!!
	return ( ( ( iT * nmub + imub ) * nmuq + imuq ) * nmus + imus );
}

bool are_nearby(const vector<double> & a, const vector<double> & b, double r0)
{
	if (mode == 1)
		return (  (a[0]-b[0])*(a[0]-b[0])
			+ (a[1]-b[1])*(a[1]-b[1])
			+ (a[2]-b[2])*(a[2]-b[2])
			+ (a[3]-b[3])*(a[3]-b[3]) < r0*r0 );
	else if (mode == 2)
		return (   (a[0]-b[0])*(a[0]-b[0]) < r0*r0
			&& (a[1]-b[1])*(a[1]-b[1]) < r0*r0
			&& (a[2]-b[2])*(a[2]-b[2]) < r0*r0
			&& (a[3]-b[3])*(a[3]-b[3]) < r0*r0 );
	else
                return (   (a[0]-b[0])*(a[0]-b[0]) < r0*r0
                        || (a[1]-b[1])*(a[1]-b[1]) < r0*r0
                        || (a[2]-b[2])*(a[2]-b[2]) < r0*r0
                        || (a[3]-b[3])*(a[3]-b[3]) < r0*r0 );

}

inline double d2( const vector<double> & a, const vector<double> & b )
{
	return ( (a[0]-b[0])*(a[0]-b[0])
			+ (a[1]-b[1])*(a[1]-b[1])
			+ (a[2]-b[2])*(a[2]-b[2])
			+ (a[3]-b[3])*(a[3]-b[3]) );
}

bool density_comp( const vector<double> & a, const vector<double> & b )
{
	if (a[0]==b[0])
	{
		if (a[1]==b[1])
		{
			if (a[2]==b[2])
			{
				if (a[3]==b[3]) return true;
				else return (a[3] < b[3]);
			}
			else return (a[2] < b[2]);
		}
		else return (a[1] < b[1]);
	}
	else return (a[0] < b[0]);
}

int main(int argc, char *argv[])
{
	// check input first
	if (argc < 2) exit(-1);

	/*cout << "Bug checks:" << endl;
	if ( point_is_in_simplex( {{0,0,0,0},{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}},
							  {0.1,0.1,0.1,0.1} ) )
		cout << " found point in this simplex!" << endl;
	else
		cout << " did not find point in this simplex!" << endl;

if (true) exit(-1);*/

	// read path to input file from command line
	string path_to_file = string(argv[1]);

	vector<double> Tvec, muBvec, muSvec, muQvec;
	vector<double> evec, bvec, svec, qvec;
	vector<vector<double> > grid;

	vector<int> Tinds(nT*nmub*nmuq*nmus), mubinds(nT*nmub*nmuq*nmus),
					muqinds(nT*nmub*nmuq*nmus), musinds(nT*nmub*nmuq*nmus);
	size_t idx = 0;
	for (int iT = 0; iT < nT; iT++)
	for (int imub = 0; imub < nmub; imub++)
	for (int imuq = 0; imuq < nmuq; imuq++)
	for (int imus = 0; imus < nmus; imus++)
	{
		Tinds[idx] = iT;
		mubinds[idx] = imub;
		muqinds[idx] = imuq;
		musinds[idx++] = imus;
	}

	// then read in file itself
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

			Tvec.push_back( Tin );
			muBvec.push_back( muBin );
			muSvec.push_back( muSin );
			muQvec.push_back( muQin );
			bvec.push_back( bin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) );		// 1/fm^3
			svec.push_back( sin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) );		// 1/fm^3
			qvec.push_back( qin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) );		// 1/fm^3
			evec.push_back( ein*Tin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) );	// MeV/fm^3

			grid.push_back( vector<double>({Tin, muBin, muQin, muSin,
											ein*Tin*Tin*Tin*Tin/(hbarc*hbarc*hbarc),
											bin*Tin*Tin*Tin/(hbarc*hbarc*hbarc),
											sin*Tin*Tin*Tin/(hbarc*hbarc*hbarc),
											qin*Tin*Tin*Tin/(hbarc*hbarc*hbarc)}) );

			if (count % 1000000 == 0) cout << "Read in " << count << " lines." << endl;
		}
	}

	infile.close();

	cout << "Initialize densities..." << endl;
	vector<vector<double> > densities;
	for ( size_t ivec = 0; ivec < evec.size(); ivec++ )
	{
		vector<double> density(5); // 5th element is cell index
		density[0] = evec[ivec];
		density[1] = bvec[ivec];
		density[2] = svec[ivec];
		density[3] = qvec[ivec];
		density[4] = ivec;
		densities.push_back( density );
	}

	double emin = *min_element(evec.begin(), evec.end());
	double emax = *max_element(evec.begin(), evec.end());
	double bmin = *min_element(bvec.begin(), bvec.end());
	double bmax = *max_element(bvec.begin(), bvec.end());
	double smin = *min_element(svec.begin(), svec.end());
	double smax = *max_element(svec.begin(), svec.end());
	double qmin = *min_element(qvec.begin(), qvec.end());
	double qmax = *max_element(qvec.begin(), qvec.end());

	// normalize to unit hypercube
	cout << "Normalizing densities..." << endl;
	for ( auto & density : densities )
	{
		density[0] = (density[0] - emin) / ( emax - emin );
		density[1] = (density[1] - bmin) / ( bmax - bmin );
		density[2] = (density[2] - smin) / ( smax - smin );
		density[3] = (density[3] - qmin) / ( qmax - qmin );
	}

	cout << "Normalizing grid cells (same thing)..." << endl;
	for ( auto & cell : grid )
	{
		cell[4] = (cell[4] - emin) / ( emax - emin );
		cell[5] = (cell[5] - bmin) / ( bmax - bmin );
		cell[6] = (cell[6] - smin) / ( smax - smin );
		cell[7] = (cell[7] - qmin) / ( qmax - qmin );
	}

	//cout << "Sorting..." << endl;
	//std::sort(densities.begin(), densities.end(), density_comp);

	const double e0 = 46308.20963821, b0 = -1.23317452, s0 = -1.53064765, q0 = -0.24540761;
	const double ne0 = (e0 - emin) / (emax - emin);
	const double nb0 = (b0 - bmin) / (bmax - bmin);
	const double ns0 = (s0 - smin) / (smax - smin);
	const double nq0 = (q0 - qmin) / (qmax - qmin);

	vector<double> nv0 = {ne0, nb0, ns0, nq0};

	Stopwatch sw;
	sw.Start();

	// get nearest neighbor w.r.t. squared distance
	vector<double> NN = *min_element(densities.begin(), densities.end(),
			[nv0](const vector<double> & a, const vector<double> & b)
			{ return d2(a, nv0) < d2(b, nv0); });

	const size_t NN_index = NN[4];
	const int iTNN = Tinds[NN_index], imubNN = mubinds[NN_index],
				imuqNN = muqinds[NN_index], imusNN = musinds[NN_index];

	sw.Stop();
	cout << "Brute force: NN_index = " << NN_index << endl;
	cout << "Brute force: NN = " << NN[0] << "   " << NN[1] << "   " << NN[2] << "   " << NN[3] << endl;
	cout << "Indentified NN simplices in " << sw.printTime() << " s." << endl;

cout << "NN grid elements are:" << endl;
	for ( auto & elem : grid[indexer( iTNN, imubNN, imuqNN, imusNN ) ] )
		cout << elem << "   ";
	cout << endl << endl;

cout << "NN grid elements are (should be the same as above):" << endl;
	for ( auto & elem : grid[NN_index] )
		cout << elem << "   ";
	cout << endl << endl;

if ( true) exit(-1);

cout << "Points in grid around NN are:" << endl;
	for (int ii = -1; ii <= 1; ii++)
        for (int jj = -1; jj <= 1; jj++)
        for (int kk = -1; kk <= 1; kk++)
        for (int ll = -1; ll <= 1; ll++)
	{
		if ( iTNN+ii < nT && iTNN+ii >= 0
		&& imubNN+jj < nmub && imubNN+jj >= 0
		&& imuqNN+kk < nmuq && imuqNN+kk >= 0
		&& imusNN+ll < nmus && imusNN+ll >= 0 )
		{
//		cout << ii << "   " << jj << "   " << kk << "   " << ll << ": ";
//		for ( auto & elem : grid[indexer( iTNN+ii, imubNN+jj, imuqNN+kk, imusNN+ll ) ] )
//			cout << elem << "   ";
//		cout << endl;
		cout << ii << "   " << jj << "   " << kk << "   " << ll << ": ";
		const auto & gridVertex = grid[indexer( iTNN+ii, imubNN+jj, imuqNN+kk, imusNN+ll ) ];
		cout << gridVertex[0] << "   " << gridVertex[1] << "   "
				<< gridVertex[2] << "   " << gridVertex[3] << "   "
				<< emin + gridVertex[4]*(emax-emin) << "   "
				<< bmin + gridVertex[5]*(bmax-bmin) << "   "
				<< smin + gridVertex[6]*(smax-smin) << "   "
				<< qmin + gridVertex[7]*(qmax-qmin) << endl;
		}
	}
	cout << endl;

//if (true) exit(-1);

	vector<vector<double> > vertices;
	for (int ii = -1; ii <= 1; ii++)
	for (int jj = -1; jj <= 1; jj++)
	for (int kk = -1; kk <= 1; kk++)
	for (int ll = -1; ll <= 1; ll++)
	{	if ( iTNN+ii < nT && iTNN+ii >= 0
			&& imubNN+jj < nmub && imubNN+jj >= 0
			&& imuqNN+kk < nmuq && imuqNN+kk >= 0
			&& imusNN+ll < nmus && imusNN+ll >= 0 )
		vertices.push_back( grid[indexer( iTNN+ii, imubNN+jj, imuqNN+kk, imusNN+ll )] );
	}

	// Qhull requires vertices as 1D vector
	vector<double> verticesFlat;
	for (int ii = -1; ii <= 1; ii++)
	for (int jj = -1; jj <= 1; jj++)
	for (int kk = -1; kk <= 1; kk++)
	for (int ll = -1; ll <= 1; ll++)
	{
		const vector<double> & gridVertex = grid[indexer( iTNN+ii, imubNN+jj, imuqNN+kk, imusNN+ll )];
		if ( iTNN+ii < nT && iTNN+ii >= 0
			&& imubNN+jj < nmub && imubNN+jj >= 0
			&& imuqNN+kk < nmuq && imuqNN+kk >= 0
			&& imusNN+ll < nmus && imusNN+ll >= 0 )
		{
		verticesFlat.push_back( gridVertex[4] );	//e
		verticesFlat.push_back( gridVertex[5] );	//b
		verticesFlat.push_back( gridVertex[6] );	//s
		verticesFlat.push_back( gridVertex[7] );	//q
		}
	}

	size_t densities_size = densities.size();
	std::vector<std::array<double, 4> > density_points(densities_size);
	for (size_t ii = 0; ii < densities_size; ii++)
		std::copy_n( densities[ii].begin(), 4, density_points[ii].begin() );

	// try this
	try
	{
		typedef point<double, 4> point4d;
		typedef kdtree<double, 4> tree4d;

		sw.Reset();
		sw.Start();
		tree4d tree(std::begin(density_points), std::end(density_points));
		sw.Stop();
		cout << "Constructed full tree in " << sw.printTime() << " s." << endl;

		sw.Reset();
		sw.Start();
		//for (size_t ii = 0; ii < 10000000; ii++)
		//	point4d n0 = tree.nearest({ne0+ii*1e-8, nb0, ns0, nq0});
		size_t kdtree_nn_index = 0;
		point4d n = tree.nearest({ne0, nb0, ns0, nq0}, kdtree_nn_index);
		sw.Stop();
		cout << "KD-Tree: Found nearest neighbor in " << setprecision(18)
				<< sw.printTime() << " s." << endl;
		cout << "KD-Tree: Nearest neighbor is " << n << endl;
		cout << "KD-Tree: Nearest neighbor index is " << kdtree_nn_index << endl;
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}



	// Test the Delaunay part here
	// first get the triangulation
	vector<vector<size_t> > simplices;
	sw.Reset();
	sw.Start();
	compute_delaunay(&verticesFlat[0], 4, verticesFlat.size() / 4, simplices);
	sw.Stop();
	cout << "Finished the Delaunay triangulation in " << sw.printTime() << " s." << endl;
	
	cout << endl << "The triangulation contains the following simplices:"  << setprecision(8) << endl;
	int isimplex = 0;
	int iclosestsimplex = 0;
	double center_d2_min = 2.0;	// start with unrealistically large value (0 <= d2 <= 1)
	for ( const auto & simplex : simplices )
	{
		cout << isimplex++ << ":";
		vector<double> center(4, 0.0);
		bool NN_vertex_included = false;
		for ( const auto & vertex : simplex )
		{
			std::transform( center.begin(), center.end(), vertices[vertex].begin()+4,
							center.begin(), std::plus<double>());
			cout << "   " << vertex;
			if ( vertex == vertices.size()/2 ) NN_vertex_included = true;
		}

		// center is average of this simplex's vertices
		std::transform( center.begin(), center.end(), center.begin(),
						[](double & element){ return 0.2*element; } );	// 0.2 == 1/(dim+1), dim == 4
		// print center
		cout << "   center: " << center[0] << "   " << center[1] << "   "
				<< center[2] << "   " << center[3];
		// print if this simplex includes nearest neighbor
		string inclusion_string = NN_vertex_included? "NN IS INCLUDED" : "NN IS NOT INCLUDED";
		cout << "; " << inclusion_string << endl;
		double d2loc = d2( center, nv0 );
		if ( d2loc < center_d2_min )
		{
			iclosestsimplex = isimplex;
			center_d2_min = d2loc;
		}
	}

	cout << "The closest simplex was isimplex = " << iclosestsimplex
			<< " with d = " << sqrt(center_d2_min) << endl;
	cout << "This simplex contains the following vertices: " << endl;
	int ivertex = 0;
	vector<vector<double> > simplexVertices;
	for ( const auto & vertex : simplices[iclosestsimplex] )
	{
		cout << ivertex++ << ":";
		for ( auto it = vertices[vertex].begin()+4; it != vertices[vertex].end(); it++ )
			cout << "   " << *it;
		cout << endl;
		simplexVertices.push_back( vector<double>( vertices[vertex].begin()+4, vertices[vertex].end() ) );
	}

	cout << endl
			<< "********************************" << endl
			<< "On to the interpolation!" << endl
			<< "********************************" << endl << endl;

	cout << "Next, check if point is contained in the simplex we found:" << endl;
	bool foundPoint = point_is_in_simplex( simplexVertices, nv0 );
	
	cout << "nv0: " << nv0[0] << "   " << nv0[1] << "   " << nv0[2] << "   " << nv0[3] << endl;


	// create some space
	for ( int iloop = 0; iloop < 10; iloop++) cout << endl;

	
	cout << "Looping through all simplices:" << endl;
	isimplex = 0;
	for ( auto & simplex : simplices )
	{
		cout << isimplex++ << ":";
		for ( const auto & vertex : simplex )
			simplexVertices.push_back( vector<double>( vertices[vertex].begin()+4,
														vertices[vertex].end() ) );
		if ( point_is_in_simplex( simplexVertices, nv0, false ) )
			cout << " found point in this simplex!" << endl;
		else
			cout << " did not find point in this simplex!" << endl;
	}


	return 0;
}

