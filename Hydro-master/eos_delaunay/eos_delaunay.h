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

		// constructors and initialization routines
		eos_delaunay(string EoS_table_file, int e_or_s);
	    eos_delaunay();
		void init(string EoS_table_file, int e_or_s);

		bool interpolate(const vector<double> & v0, vector<double> & result, bool verbose = false);
		bool interpolate_NNmode(const vector<double> & v0, vector<double> & result, bool verbose = false);
		bool interpolate_NMNmode(const vector<double> & v0, vector<double> & result, bool verbose = false);
		bool interpolate_NMNmode_v2(const vector<double> & v0, vector<double> & result,
									bool expand_hypercube = false, bool verbose = false);
		bool interpolate_NMNmode_v3(const vector<double> & v0, vector<double> & result, bool verbose = false);
		bool triangulate_and_locate_point( const vector<double> & nv0, const vector<int> & base,
					vector<vector<double> > & vertices, vector<vector<size_t> > & simplices,
					vector<double> & point_lambda_in_simplex, int & iclosestsimplex );
		void refine_hypercube(vector<vector<double> > & hypercube);

	private:

		typedef point<double, 4> point4d;
		typedef kdtree<double, 4> tree4d;

		tree4d * tree_ptr, * midpoint_tree_ptr, * unnormalized_midpoint_tree_ptr;

		static const double hbarc;
		//const size_t nT = 155, nmub = 37, nmus = 37, nmuq = 37;
		static const size_t nT, nmub, nmus, nmuq;

		double emin, emax, bmin, bmax, smin, smax, qmin, qmax;

		int using_e_or_s_mode;

		vector<vector<double> > grid, unnormalized_grid;
		vector<int> Tinds, mubinds, muqinds, musinds;
		vector<vector<size_t> > midpoint_inds;
		
		inline size_t indexer( const int iT, const int imub, const int imuq, const int imus )
		{
			// mus varies faster than muq!!!!!
			return ( ( ( iT * nmub + imub ) * nmuq + imuq ) * nmus + imus );
		}
//		size_t indexer( const int iT, const int imub, const int imuq, const int imus );

		
		inline double d2( const vector<double> & a, const vector<double> & b )
		{
			return (  (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1])
					+ (a[2]-b[2])*(a[2]-b[2]) + (a[3]-b[3])*(a[3]-b[3]) );
		}
		
		
		// prototypes
		void load_EoS_table(string path_to_file, vector<vector<double> > & grid, int e_or_s);
		void get_min_and_max(vector<double> & v, double & minval, double & maxval, bool normalize);

};



#endif
