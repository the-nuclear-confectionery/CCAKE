#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string>

#include <lib.h>
#include "eos_delaunay.h"
#include "Stopwatch.h"

using namespace std;

void load_EoS_table(string path_to_file, vector<vector<double> > & grid);

int main(int argc, char *argv[])
{
	// check input first
	if (argc < 3) exit(-1);

	constexpr bool timing_test_only = false;
	Stopwatch sw;

	/*if (!timing_test_only)
	{
		cout << "Testing interface to C library containing BSQ EoS" << endl;
		printf("Starting!\n");
	        initialize("../Coefficients_Parameters.dat");
	
	        double point[4] = {200.0, 50.0, 75.0, 100.0};
	        double densities[4];
	        get_densities(point, densities);
	        printf("Found these densities: %lf %lf %lf %lf\n",
	                densities[0], densities[1], densities[2], densities[3]);
	
	        printf("Finished!\n");
		//if (true) exit(-1);
	}*/

	sw.Reset();
	sw.Start();

	// read path to input file from command line
	string path_to_file = string(argv[1]);
	string path_to_staggered_file = string(argv[2]);

	// set up EoS object
	initialize("../Coefficients_Parameters.dat");	// this is the C library that must be initialized be using
	eos_delaunay EoS( path_to_file );

	size_t cellCount = 0;

	vector<double> result(4, 0.0);
	sw.Reset();
	sw.Start();
	//EoS.interpolate({3405.08, -0.473819, -1.78269, -2.89511}, result);
	//EoS.interpolate({1940.68, -0.284676, -1.0705, -1.5329}, result);
	//EoS.interpolate({1152.45, -0.402379, 0.562341, 0.00269458}, result);	// Exact: 152.5   -387.5   -37.5   212.5
	//EoS.interpolate({973.563, -0.316059, 0.323859, 1.06384}, result);
	//EoS.interpolate({31.0278, -8.98404e-05, -2.03722e-05, -0.0437323}, result);
	EoS.interpolate({573.228, -0.0236158, -0.129289, -0.711049}, result, true);
	cout << "Exact: " << "102.5   -437.5   -437.5   -212.5" << endl;
	sw.Stop();
	cout << "Found the solution in " << sw.printTime() << " s." << endl;
	for (const double & elem : result)
				cout << "   " << elem;
			cout << endl;
	
	if (true) exit(-1);
	
	// load staggered file with test points
	vector<vector<double> > staggered_grid;
	load_EoS_table(path_to_staggered_file, staggered_grid);


	sw.Stop();
	cout << "Set-up took " << sw.printTime() << " s." << endl;

	sw.Reset();
	sw.Start();
	for (const vector<double> & staggered_cell : staggered_grid)
	{
		if ( staggered_cell[0] < 100.0 || staggered_cell[0] > 200.0 )	// MeV
			continue;
		// interpolate the densities
		EoS.interpolate(vector<double>(staggered_cell.begin()+4,staggered_cell.end()), result);
		if ( timing_test_only )
		{
			if (++cellCount % 1000 == 0)
			{
				sw.Stop();
				cout << "Interpolated " << cellCount << " cells in " << sw.printTime() << " s." << endl;
				//sw.Reset();
				sw.Start();
			}
		}
		else
		{
			cout << cellCount++ << ":";
			for (const double & elem : staggered_cell)
				cout << "   " << elem;
			for (const double & elem : result)
				cout << "   " << elem;
			cout << "\n";
			//if (true){ cout << endl; exit(-1); }
		}
	}

	/*
	// set densities where we want to test the interpolator
	const double e0 = 46308.20963821, b0 = -1.23317452, s0 = -1.53064765, q0 = -0.24540761;

	// Time set-up
	Stopwatch sw;
	sw.Reset();
	sw.Start();

	vector<double> result(4);
	eos_delaunay EoS( path_to_file );

	sw.Stop();
	cout << "Set-up took " << sw.printTime() << " s." << endl;

	sw.Reset();
	sw.Start();

	EoS.interpolate({e0, b0, s0, q0}, result);

	sw.Stop();

	cout << "The final answer is: "
		<< result[0] << "   " << result[1] << "   " << result[2] << "   " << result[3] << endl;
		//<< T0 << "   " << mub0 << "   " << muq0 << "   " << mus0 << endl;
	cout << "Found the answer in approximately " << sw.printTime() << " s." << endl;
	*/

	return 0;
}

void load_EoS_table(string path_to_file, vector<vector<double> > & grid)
{
	const double hbarc = 197.327;
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

			grid.push_back( vector<double>({Tin, muBin, muQin, muSin,
											ein*Tin*Tin*Tin*Tin/(hbarc*hbarc*hbarc),
											bin*Tin*Tin*Tin/(hbarc*hbarc*hbarc),
											sin*Tin*Tin*Tin/(hbarc*hbarc*hbarc),
											qin*Tin*Tin*Tin/(hbarc*hbarc*hbarc)}) );

			if (count % 1000000 == 0)
				cout << "Read in " << count << " lines of "
					<< path_to_file << "." << endl;
		}
	}

	infile.close();
	return;
}

