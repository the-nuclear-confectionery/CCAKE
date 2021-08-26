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

inline size_t indexer( const int iT, const int imub, const int imuq, const int imus )
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


// prototypes
void load_EoS_table(std::string path_to_file, vector<vector<double> > & grid);
void get_min_and_max(vector<double> & v, double & minval, double & maxval, bool normalize);


int main(int argc, char *argv[])
{
	// check input first
	if (argc < 2) exit(-1);

	// read path to input file from command line
	string path_to_file = string(argv[1]);

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
		musinds[idx] = imus;
		idx++;
	}

	// load EoS table
	load_EoS_table(path_to_file, grid);

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
	double emin = 0.0, emax = 0.0, bmin = 0.0, bmax = 0.0,
			smin = 0.0, smax = 0.0, qmin = 0.0, qmax = 0.0;
	get_min_and_max(evec, emin, emax, false);
	get_min_and_max(bvec, bmin, bmax, false);
	get_min_and_max(svec, smin, smax, false);
	get_min_and_max(qvec, qmin, qmax, false);

	// normalize grid points
	for ( vector<double> & gridcell : grid )
	{
		gridcell[4] = (gridcell[4] - emin) / ( emax - emin );
		gridcell[5] = (gridcell[5] - bmin) / ( bmax - bmin );
		gridcell[6] = (gridcell[6] - smin) / ( smax - smin );
		gridcell[7] = (gridcell[7] - qmin) / ( qmax - qmin );
	}

	// copy normalized densities from grid to separate vector (for kd-tree)
	// (needs to be vector of arrays to set up kd-tree correctly)
	std::vector<std::array<double, 4> > density_points(grid.size());
	for (size_t ii = 0; ii < grid.size(); ii++)
		std::copy_n( grid[ii].begin()+4, 4, density_points[ii].begin() );

	// set up kd-tree
	typedef point<double, 4> point4d;
	typedef kdtree<double, 4> tree4d;
	//try
	//{
		cout << "Setting up kd-tree...";
		tree4d tree(std::begin(density_points), std::end(density_points));
		cout << "finished!\n";
		//cout << "Constructed full tree in " << sw.printTime() << " s." << endl;
	//}
	//catch (const std::exception& e)
	//{
	//	std::cerr << '\n' << e.what() << '\n';
	//}

	// set densities where we want to test the interpolator
	//const double e0 = 46308.20963821, b0 = -1.23317452, s0 = -1.53064765, q0 = -0.24540761;
	const double e0 = 3405.08, b0 = -0.473819, s0 = -1.78269, q0 = -2.89511;

	// Set-up is finished; start timing now
	Stopwatch sw;
	sw.Reset();
	sw.Start();

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
		point4d n = tree.nearest({ne0, nb0, ns0, nq0}, kdtree_nn_index);
		//sw.Stop();
		cout << "Query point: {" << ne0 << ", " << nb0 << ", " << ns0 << ", " << nq0 << "}" << endl;
		cout << "KD-Tree: Found nearest neighbor in " << setprecision(18)
				<< sw.printTime() << " s." << endl;
		cout << "KD-Tree: Nearest neighbor is " << n << endl;
		cout << "KD-Tree: Nearest neighbor index is " << kdtree_nn_index << endl;
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
	compute_delaunay(&verticesFlat[0], 4, verticesFlat.size() / 4, simplices);
	//cout << "Finished the Delaunay triangulation in " << sw.printTime() << " s." << endl;


	// =======================================================
	// triangulation is complete; now find containing simplex

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


	cout << "simplices[iclosestsimplex].size() = " << simplices[iclosestsimplex].size() << endl;
	{
		int ivertex = 0;
	for ( const auto & vertex : simplexVertices )
	{
		cout << "ivertex = " << ivertex++ << ":" << endl;
		for ( const auto & coordinate : vertex )
			cout << "   " << coordinate;
		cout << endl;
	}
	}

	// locate the point in the simplex (assuming we know it's there;
	// currently no plan B for if point winds up outside this simplex)
	vector<double> point_lambda_in_simplex(5, 0.0);	// dim + 1 == 5
	bool foundPoint = point_is_in_simplex( simplexVertices, nv0, point_lambda_in_simplex, true );


	// finally, use the output lambda coefficients to get the interpolated values
	double T0 = 0.0, mub0 = 0.0, muq0 = 0.0, mus0 = 0.0;
	{
		int ivertex = 0;
		for ( const auto & vertex : simplices[iclosestsimplex] )
		{
			double lambda_coefficient = point_lambda_in_simplex[ivertex];
			cout << "Check interpolation: " << ivertex << "   " << vertex << "   "
					<< point_lambda_in_simplex[ivertex] << "   "
					<< vertices[vertex][0] << "   " << vertices[vertex][1] << "   " 
					<< vertices[vertex][2] << "   " << vertices[vertex][3] << endl;
			T0   += lambda_coefficient * vertices[vertex][0];
			mub0 += lambda_coefficient * vertices[vertex][1];
			muq0 += lambda_coefficient * vertices[vertex][2];
			mus0 += lambda_coefficient * vertices[vertex][3];
			ivertex++;
		}
	}

	sw.Stop();

	cout << "The final answer is: "
		<< T0 << "   " << mub0 << "   " << muq0 << "   " << mus0 << endl;

	cout << "Found the answer in approximately " << sw.printTime() << " s." << endl;

	return 0;
}




void load_EoS_table(std::string path_to_file, vector<vector<double> > & grid)
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

			if (count % 1000000 == 0) cout << "Read in " << count << " lines." << endl;
		}
	}

	infile.close();
	return;
}

void get_min_and_max(vector<double> & v, double & minval, double & maxval, bool normalize)
{
	minval = *min_element(v.begin(), v.end());
	maxval = *max_element(v.begin(), v.end());
	if (normalize)
		std::transform( v.begin(), v.end(), v.begin(),
						[minval,maxval](double & element)
						{ return (element-minval)/(maxval-minval); } );
	return;
}
