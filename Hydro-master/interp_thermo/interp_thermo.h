#ifndef INTERP_THERMO_H
#define INTERP_THERMO_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <gsl/gsl_linalg.h>

#include "Stopwatch.h"

namespace interp_thermo
{
	double emin, emax, rhoBmin, rhoBmax, rhoSmin, rhoSmax, rhoQmin, rhoQmax;
	
	// probably don't need a separate function to do this...
	inline double normalize( double x0, double x1, double x )
	{
		return ( (x-x0)/(x1-x0) );
	}
	
	// distances computed in density space (last 4 of 8 elements in first vector)
	// ==>> fixed point must be second argument
	inline double distance2( const vector<double> & a, const vector<double> & b )
	{
		return (a[4]-b[0])*(a[4]-b[0]) + (a[5]-b[1])*(a[5]-b[1])
			 + (a[6]-b[2])*(a[6]-b[2]) + (a[7]-b[3])*(a[7]-b[3]);
	}
	
	// read in EoS table to be interpolated over and normalize densities to respective ranges
	void load_file( string filename, vector<vector<double> > & EoS_table )
	{
		Stopwatch sw;
		sw.Start();
		const double hbarc = 197.33;
		const double hbarc3 = hbarc*hbarc*hbarc;
		EoS_table.clear();
		ifstream infile(filename.c_str());
		if (infile.is_open())
		{
			string line;
			int count = 0;
			double T, muB, muQ, muS, e, rhoB, rhoS, rhoQ, dummy;
			while ( getline (infile, line) )
			{
				istringstream iss(line);
				iss >> T >> muB >> muQ >> muS
					>> dummy >> dummy
					>> rhoB >> rhoS >> rhoQ >> e;
	
				e    *= T*T*T*T/hbarc3;	// MeV/fm^3
				rhoB *= T*T*T/hbarc3;	// 1/fm^3
				rhoS *= T*T*T/hbarc3;	// 1/fm^3
				rhoQ *= T*T*T/hbarc3;	// 1/fm^3
	
				vector<double> EoS_entry;
				EoS_entry.push_back(T);
				EoS_entry.push_back(muB);
				EoS_entry.push_back(muS);
				EoS_entry.push_back(muQ);
				EoS_entry.push_back(e);
				EoS_entry.push_back(rhoB);
				EoS_entry.push_back(rhoS);
				EoS_entry.push_back(rhoQ);
	
				EoS_table.push_back( EoS_entry );
	
				if ( count++ < 1 )
				{
					emin = e; emax = e;
					rhoBmin = rhoB; rhoBmax = rhoB;
					rhoSmin = rhoS; rhoSmax = rhoS;
					rhoQmin = rhoQ; rhoQmax = rhoQ;
				}
				else
				{
					if ( e < emin ) emin = e;
					if ( e > emax ) emax = e;
					if ( rhoB < rhoBmin ) rhoBmin = rhoB;
					if ( rhoB > rhoBmax ) rhoBmax = rhoB;
					if ( rhoS < rhoSmin ) rhoSmin = rhoS;
					if ( rhoS > rhoSmax ) rhoSmax = rhoS;
					if ( rhoQ < rhoQmin ) rhoQmin = rhoQ;
					if ( rhoQ > rhoQmax ) rhoQmax = rhoQ;
				}
			}
		}
	
		infile.close();
	
		// normalize input data (easy to go the other direction too)
		for ( auto & EoS_entry : EoS_table )
		{
			EoS_entry[4] = normalize( emin,    emax,    EoS_entry[4] );
			EoS_entry[5] = normalize( rhoBmin, rhoBmax, EoS_entry[5] );
			EoS_entry[6] = normalize( rhoSmin, rhoSmax, EoS_entry[6] );
			EoS_entry[7] = normalize( rhoQmin, rhoQmax, EoS_entry[7] );
		}
	
		sw.Stop();
		cout << "Finished loading " << filename << " in " << sw.printTime() << " s." << endl;
	
		return;
	}
	
	// compute nearest neighbors with a partial sort and copy entries to vector
	void get_nearest_neighbors( const vector<vector<double> > & EoS_table,
								vector<vector<double> > & neighbors,
								const vector<double> & p, const size_t k )
	{
		neighbors.clear();
		neighbors = vector<vector<double> > ( k, vector<double> ( p.size(), 0.0 ) );
	
		partial_sort_copy(
			EoS_table.begin(), EoS_table.end(), neighbors.begin(), neighbors.end(),
			[&p](const vector<double> & x, const vector<double> & y) -> bool
				{ return distance2(x,p) < distance2(y,p); } );
	
		return;
	}
	
	// given vertices, define simplex and get barycentric w.r.t. it for additional point
	void get_barycentric_coordinates(double T[], double d[], double lambda[], int dim)
	{
		gsl_matrix_view m = gsl_matrix_view_array (T, dim, dim);
		gsl_vector_view b = gsl_vector_view_array (d, dim);
		gsl_vector *x = gsl_vector_alloc (dim);
		
		int s;
		gsl_permutation * p = gsl_permutation_alloc (dim);
		
		gsl_linalg_LU_decomp (&m.matrix, p, &s);
		gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
	
		for (int i = 0; i < dim; i++) lambda[i] = gsl_vector_get(x, i);
		
		gsl_permutation_free (p);
		gsl_vector_free (x);
	
		return;
	}
}

#endif
