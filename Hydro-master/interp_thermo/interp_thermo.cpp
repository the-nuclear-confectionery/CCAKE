#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "interp_thermo.h"

using namespace std;

// prototypes
//void get_barycentric_coordinates(double T[], double d[], double lambda[], int dim);
//void load_file( string filename, vector<vector<double> > & EoS_table );
//void get_nearest_neighbors( const vector<vector<double> > & EoS_table,
//							vector<vector<double> > & neighbors,
//							const vector<double> & p, const size_t k );

// driver function
int main(int argc, char ** argv)
{
	// load thermo file
	string filename = argv[1];
	vector<vector<double> > EoS_table;
	interp_thermo::load_file(filename, EoS_table);

	// pick a point
	vector<double> point {0.5, 0.5, 0.5, 0.5};  // normalized coordinates

	// get k nearest neighbors to point
	const size_t k = 10;
	vector<vector<double> > neighbors;
	interp_thermo::get_nearest_neighbors( EoS_table, neighbors, point, k );

	cout << "Nearest neighbors are:" << endl;
	for ( auto & neighbor : neighbors )
	{
		for ( auto & element : neighbor )
			cout << element << "   ";
		cout << endl;
	}
	
	return 0;
}

