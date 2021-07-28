#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

double emin, emax, rhoBmin, rhoBmax, rhoSmin, rhoSmax, rhoQmin, rhoQmax;

// prototypes
void load_file( string filename, vector<vector<double> > & EoS_table );
void get_nearest_neighbors( const vector<vector<double> > & EoS_table,
							vector<vector<double> > & neighbors,
							const vector<double> & p, const int k );

// driver function
int main(int argc, char ** argv)
{
	// load thermo file
	string filename = argv[1];
	vector<vector<double> > EoS_table;
	load_file(filename, EoS_table);

	// pick a point
	vector<double> point {0.5, 0.5, 0.5, 0.5};  // normalized coordinates

	// get k nearest neighbors to point
	const int k = 10;
	vector<vector<double> > neighbors;
	get_nearest_neighbors( EoS_table, neighbors, point, k );

	cout << "Nearest neighbors are:" << endl;
	for ( auto & neighbor : neighbors )
	{
		for ( auto & element : neighbor )
			cout << element << "   ";
		cout << endl;
	}
	
	return 0;
}

inline double normalize( double x0, double x1, double x )
{
	return ( (x-x0)/(x1-x0) );
}

void load_file( string filename, vector<vector<double> > & EoS_table )
{
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

	return;
}

// distances computed in density space (last 4 of 8 elements in vector)
inline double distance2( const vector<double> & a, const vector<double> & b )
{
	return (a[4]-b[4])*(a[4]-b[4]) + (a[5]-b[5])*(a[5]-b[5])
		 + (a[6]-b[6])*(a[6]-b[6]) + (a[7]-b[7])*(a[7]-b[7]);
}

void get_nearest_neighbors( const vector<vector<double> > & EoS_table,
							vector<vector<double> > & neighbors,
							const vector<double> & p, const int k )
{
	neighbors.clear();
	neighbors = vector<vector<double> > ( k, vector<double> ( p.size(), 0.0 ) );

	bool closer_than = [&p](const vector<double> & x, const vector<double> & y)
					   { return distance2(x,p) < distance2(y,p); };

	partial_sort_copy( EoS_table.begin(), EoS_table.end(),
						neighbors.begin(), neighbors.end(), closer_than );

	return;
}

