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
	constexpr int test_mode = 2;

	if ( test_mode == 1 )	// check single point
	{
		// load thermo file
		string filename = argv[1];
		vector<vector<double> > EoS_table;
		interp_thermo::load_file(filename, EoS_table);
	
		// pick a point
		vector<double> point {0.5, 0.5, 0.5, 0.5};  // normalized coordinates
	
		// get k nearest neighbors to point
		const size_t k = 81;
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
	}
	else	// read in all points to test from file
	{
		// load thermo file
		string filename = argv[1];
		string testfilename = argv[2];
		vector<vector<double> > EoS_table, points_to_check;
		interp_thermo::load_file(filename, EoS_table);
		interp_thermo::load_test_file(testfilename, points_to_check);

		int count = 0;
		size_t k = 4;
		do
		{
			cout << "k = " << k << endl;
			vector<vector<double> > neighbors;
			for ( auto & point_to_check : points_to_check )
			{

			cout << "Point is at:" << endl;
			for ( auto & element : point_to_check )
					cout << element << "   ";
				cout << endl;

				interp_thermo::get_nearest_neighbors( EoS_table, neighbors, point_to_check, k, true );
	
			cout << "Nearest neighbors are:" << endl;
			for ( auto & neighbor : neighbors )
			{
				for ( auto & element : neighbor )
					cout << element << "   ";
				cout << endl;
			}

				for ( double power_in = 1.0; power_in <= 5.0; power_in += 0.5 )
					interp_thermo::get_IDW_point_estimate( neighbors, point_to_check, power_in );
			}
			k *= 2;
			cout << endl;
		} while ( count++ < 1 );
		
	}
	
	return 0;
}

