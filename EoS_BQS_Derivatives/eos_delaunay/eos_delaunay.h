#ifndef EOS_DELAUNAY_H
#define EOS_DELAUNAY_H

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

using namespace std;

class eos_delaunay
{
	public:

		eos_delaunay(string EoS_table_file);
		//eos_delaunay(vector<vector<double> > & grid_in);

	    //eos_delaunay();

	    //void init(string EoS_table_file);
	    //void init(vector<vector<double> > & grid_in);

		void interpolate(const vector<double> & v0, vector<double> & result);
		bool interpolate_NMNmode(const vector<double> & v0, vector<double> & result);
		bool interpolate_NMNmode_v2(const vector<double> & v0, vector<double> & result);
		bool triangulate_and_locate_point( const vector<double> & nv0, const vector<int> & base,
					vector<vector<double> > & vertices, vector<vector<size_t> > & simplices,
					vector<double> & point_lambda_in_simplex, int & iclosestsimplex );
		void refine_hypercube(vector<vector<double> > & hypercube);

	private:

		typedef point<double, 4> point4d;
		typedef kdtree<double, 4> tree4d;

		tree4d * tree_ptr, * midpoint_tree_ptr, * unnormalized_midpoint_tree_ptr;

		const double hbarc = 197.327;
		const size_t nT = 155, nmub = 37, nmus = 37, nmuq = 37;

		double emin, emax, bmin, bmax, smin, smax, qmin, qmax;

		vector<vector<double> > grid, unnormalized_grid;
		vector<int> Tinds, mubinds, muqinds, musinds;
		vector<vector<size_t> > midpoint_inds;
		
		inline size_t indexer( const int iT, const int imub, const int imuq, const int imus )
		{
			// mus varies faster than muq!!!!!
			return ( ( ( iT * nmub + imub ) * nmuq + imuq ) * nmus + imus );
		}		
		
		/*inline size_t indexer( const vector<int> & inds )
		{
			// mus varies faster than muq!!!!!
			return ( ( ( inds[0] * nmub + inds[1] ) * nmuq + inds[2] ) * nmus + inds[3] );
		}*/	
		
		inline double d2( const vector<double> & a, const vector<double> & b )
		{
			return (  (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1])
					+ (a[2]-b[2])*(a[2]-b[2]) + (a[3]-b[3])*(a[3]-b[3]) );
		}
		
		
		// prototypes
		void load_EoS_table(string path_to_file, vector<vector<double> > & grid);
		void get_min_and_max(vector<double> & v, double & minval, double & maxval, bool normalize);

};



#endif
