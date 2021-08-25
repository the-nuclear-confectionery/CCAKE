#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string>

//#include "delaunay.h"
//#include "kdtree.h"
//#include "point_in_simplex.h"
#include "eos_delaunay.h"
#include "Stopwatch.h"

using namespace std;

int main(int argc, char *argv[])
{
	// check input first
	if (argc < 2) exit(-1);

	// read path to input file from command line
	string path_to_file = string(argv[1]);

	// set densities where we want to test the interpolator
	const double e0 = 46308.20963821, b0 = -1.23317452, s0 = -1.53064765, q0 = -0.24540761;

	// Set-up is finished; start timing now
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

	return 0;
}

