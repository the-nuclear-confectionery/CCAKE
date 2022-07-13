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
void load_ICCING(string path_to_file, vector<vector<double> > & points);

int main(int argc, char *argv[])
{
	// check input first
	if (argc < 4)
	{
		cerr << "Usage: ./interpolate_ebsq [filename of EoS table] [filename to read in] [execution mode]\n";
		exit(-1);
	}

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
	int execution_mode = std::stoi(argv[3]);

	// set up EoS object
	initialize("../Coefficients_Parameters.dat");	// this is the C library that must be initialized be using
	eos_delaunay EoS( path_to_file );

	size_t cellCount = 0;

	vector<double> result(4, 0.0);
	if ( execution_mode == 0 )
	{
	sw.Reset();
	sw.Start();
	//EoS.interpolate({3405.08, -0.473819, -1.78269, -2.89511}, result);
	//EoS.interpolate({1940.68, -0.284676, -1.0705, -1.5329}, result);
	//EoS.interpolate({1152.45, -0.402379, 0.562341, 0.00269458}, result);	// Exact: 152.5   -387.5   -37.5   212.5
	//EoS.interpolate({973.563, -0.316059, 0.323859, 1.06384}, result);
	//EoS.interpolate({31.0278, -8.98404e-05, -2.03722e-05, -0.0437323}, result);
	//EoS.interpolate({573.228, -0.0236158, -0.129289, -0.711049}, result, true);
	//cout << "Exact: " << "102.5   -437.5   -437.5   -212.5" << endl;
	//EoS.interpolate({5754.35, 0.00231029, 0.351709, 0.378919}, result, true);
	//cout << "Exact: " << "252.5   52.5   52.5   52.5" << endl;
	EoS.interpolate({75633.1, 0.0078043, -0.00259725, 0.0181128}, result, true);
	/*double e0 = 75633.1, b0 = 0.0078043, s0 = -0.00259725, q0 = 0.0181128;
	double taufactor = 1.0, Delta_tau = 1.0, taumax = 1000.0;
	while ( result[0] < 1e-6 && taufactor <= taumax ) // no temperatures this small
	{
		taufactor += Delta_tau;
		EoS.interpolate({e0/taufactor, b0/taufactor, s0/taufactor, q0/taufactor}, result, false);
	
	for (const double & elem : result)
				cout << "   " << elem;
			cout << endl;

	cout << "taufactor = " << taufactor << endl << endl;
	}*/

        for (const double & elem : result)
                                cout << "   " << elem;
                        cout << endl;
	
	sw.Stop();
	cout << "Found the solution in " << sw.printTime() << " s." << endl;
	if (true) exit(-1);
	}

	if ( execution_mode == 1 )
	{
	string path_to_staggered_file = string(argv[2]);

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
	}

	if ( execution_mode == 2 )
	{
		cout << "Loading ICCING cells\n";

		// load staggered file with test points
		vector<vector<double> > ICCING_points;
		string path_to_ICCING_file = string(argv[2]);
		load_ICCING(path_to_ICCING_file, ICCING_points);
	
	
		sw.Stop();
		cout << "Set-up took " << sw.printTime() << " s." << endl;
	
		sw.Reset();
		sw.Start();
		for (const vector<double> & ICCING_point : ICCING_points)
		{
			cellCount++;
			//if (ICCING_point[2] > 1500.0) continue;	// trying just superdense EoS

			// interpolate the densities
			EoS.interpolate(vector<double>(ICCING_point.begin()+2,ICCING_point.end()), result);

			// check solution
			double densities_arr[4];
			iter_swap(result.begin() + 2, result.begin() + 3);      // fix this eventually
			get_eBSQ_densities(result.data(), densities_arr);
			iter_swap(result.begin() + 2, result.begin() + 3);      // fix this eventually
		        vector<double> densities(densities_arr, densities_arr+4);

			cout << cellCount-1 << ":";
			for (const double & elem : ICCING_point)
				cout << "   " << elem;
			for (const double & elem : result)
				cout << "   " << elem;
			for (const double & elem : densities)
				cout << "   " << elem;
			cout << "\n";
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


void load_ICCING(string path_to_file, vector<vector<double> > & points)
{
	points.clear();
	const int n_header = 1;	// assume number of header lines is 1

	// read in file itself
	ifstream infile(path_to_file.c_str());
	if (infile.is_open())
	{
		size_t count = 0;
		string line;
		double xpt, ypt, bin, sin, qin, ein; 
		while ( getline (infile, line) )
		{
			count++;
			if (count < n_header+1) continue;

			istringstream iss(line);
			iss >> xpt >> ypt >> ein >> bin >> sin >> qin;

			points.push_back( vector<double>({xpt, ypt, 1000.0*ein, bin, sin, qin}) );

		}
	}

	infile.close();
}
