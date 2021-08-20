#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string>

#include "kdtree.h"
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

			grid.push_back( vector<double>({Tin, muBin, muQin, muSin, ein, bin, sin, qin}) );

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
	cout << "Normalizing..." << endl;
	for ( auto & density : densities )
	{
		density[0] = (density[0] - emin) / ( emax - emin );
		density[1] = (density[1] - bmin) / ( bmax - bmin );
		density[2] = (density[2] - smin) / ( smax - smin );
		density[3] = (density[3] - qmin) / ( qmax - qmin );
	}

	cout << "Sorting..." << endl;
	std::sort(densities.begin(), densities.end(), density_comp);

	Stopwatch sw;
	sw.Start();

	const double e0 = 46308.20963821, b0 = -1.23317452, s0 = -1.53064765, q0 = -0.24540761;
	const double ne0 = (e0 - emin) / (emax - emin);
	const double nb0 = (b0 - bmin) / (bmax - bmin);
	const double ns0 = (s0 - smin) / (smax - smin);
	const double nq0 = (q0 - qmin) / (qmax - qmin);

	vector<double> nv0 = {ne0, nb0, ns0, nq0};

	// get nearest neighbor w.r.t. squared distance
	vector<double> NN = *min_element(densities.begin(), densities.end(),
			[nv0](const vector<double> & a, const vector<double> & b)
			{ return d2(a, nv0) < d2(b, nv0); });

	const size_t NN_index = NN[4];
	const int iTNN = Tinds[NN_index], imubNN = mubinds[NN_index],
				imuqNN = muqinds[NN_index], imusNN = musinds[NN_index];

	vector<vector<double> > vertices;
	for (int ii = -1; ii <= 1; ii++)
	for (int jj = -1; jj <= 1; jj++)
	for (int kk = -1; kk <= 1; kk++)
	for (int ll = -1; ll <= 1; ll++)
		vertices.push_back( grid[indexer( iTNN+ii, imubNN+jj, imuqNN+kk, imusNN+ll )] );

	sw.Stop();
	cout << "Brute force: NN_index = " << NN_index << endl;
	cout << "Indentified NN simplices in " << sw.printTime() << " s." << endl;

	std::vector<std::array<double, 4>> density_points(densities_size);
	for (size_t ii = 0; ii < densities_size; ii++)
		std::copy_n( densities[ii].begin(), 4, density_points[ii].begin() );

	// try this
	try
	{
		typedef point<int, 2> point2d;
		typedef kdtree<int, 2> tree2d;
		
		point2d points[] = { { 2, 3 }, { 5, 4 }, { 9, 6 }, { 4, 7 }, { 8, 1 }, { 7, 2 } };
		
		tree2d tree(std::begin(points), std::end(points));
		point2d n = tree.nearest({ 9, 2 });
		
		std::cout << "Wikipedia example data:\n";
		std::cout << "nearest point: " << n << '\n';
		std::cout << "distance: " << tree.distance() << '\n';
		std::cout << "nodes visited: " << tree.visited() << '\n';

	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}


	return 0;
}

