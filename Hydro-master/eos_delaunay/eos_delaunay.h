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

#include "kdtree.h"

using namespace std;

class eos_delaunay
{
	public:

		// constructors and initialization routines
		eos_delaunay(string EoS_table_file, int e_or_s);
	    eos_delaunay();
		void init(string EoS_table_file, int e_or_s);

		void get_NMN_coordinates(const vector<double> & v0, vector<double> & result,
								 bool use_normalized);

		/*double get_emax() {return emax;}
		double get_emin() {return emin;}
		double get_bmax() {return bmax;}
		double get_bmin() {return bmin;}
		double get_smax() {return smax;}
		double get_smin() {return smin;}
		double get_qmax() {return qmax;}
		double get_qmin() {return qmin;}*/

		double normalized_d2(double a[], double b[]);
		double unnormalized_d2(double a[], double b[]);
		double log_d2(double a[], double b[]);

	private:

		typedef point<double, 4> point4d;
		typedef kdtree<double, 4> tree4d;

		tree4d * tree_ptr, * midpoint_tree_ptr, * unnormalized_tree_ptr,
			   * unnormalized_midpoint_tree_ptr;
		tree4d * e_tree_ptr, * e_midpoint_tree_ptr, * e_unnormalized_tree_ptr,
			   * e_unnormalized_midpoint_tree_ptr;
		tree4d * entr_tree_ptr, * entr_midpoint_tree_ptr, * entr_unnormalized_tree_ptr,
			   * entr_unnormalized_midpoint_tree_ptr;

		static const double hbarc;
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
		inline size_t indexer( const vector<size_t> & v )
		{
			// mus varies faster than muq!!!!!
			return ( ( ( v[0] * nmub + v[1] ) * nmuq + v[2] ) * nmus + v[3] );
		}
		
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
