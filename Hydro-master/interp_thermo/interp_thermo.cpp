#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "interp_thermo.h"

using namespace std;

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
	const size_t k = 25;
	vector<vector<double> > neighbors;
	interp_thermo::get_nearest_neighbors( EoS_table, neighbors, point, k );

	cout << "Nearest neighbors are:" << endl;
	for ( auto & neighbor : neighbors )
	{
		for ( auto & element : neighbor )
			cout << element << "   ";
		cout << endl;
	}

	// IGNORE THIS VERSION FOR NOW
	// now check if point is inside simplex defined by neighbors
	//interp_thermo::check_point_in_simplex( neighbors, point );

	// USE INVERSE DISTANCE WEIGHTING INSTEAD
	interp_thermo::get_IDW_point_estimate( neighbors, point );
	
	return 0;
}

