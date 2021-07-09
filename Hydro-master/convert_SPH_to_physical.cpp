#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

constexpr double PI 	= 3.1415926535897932384626433;
constexpr double h 		= 0.3;
constexpr double knorm 	= 10.0/(7.*PI*h*h);

double kernel( double x, double y )
{
	double r = sqrt(x*x+y*y);
	double q = r/h;

	if(q>=2)
		return 0;
	else if(q>=1)
		return 0.25*knorm*(2.0-q)*(2.0-q)*(2.0-q);
	else
	{
		double qq = q*q;
		return knorm*(1.0 - 1.5*qq + 0.75*q*qq);
	}
}


int main( int argc, char ** argv )
{
	string resultsDirectory 	= "./outputfiles";
	string fileStemToReadIn 	= argv[1];
	string fileStemToPrintTo 	= argv[2];
	string indexToProcess 		= argv[3];

	string infilename = resultsDirectory + "/" + fileStemToReadIn + indexToProcess + "_ev0.dat";
	string outfilename = resultsDirectory + "/" + fileStemToPrintTo + indexToProcess + "_ev0.dat";

	vector<double> xvec, yvec, pvec, Tvec, muBvec, muQvec, muSvec, evec, Bvec, Svec, Qvec, svec;
	double tau, x, y, pin, Tin, muBin, muSin, muQin, ein, Bin, Sin, Qin, sin;

	// read in data here
	int linecount = 0;
	long long nSPH = 0;
	double dummy = 0.0;
	ifstream infile( infilename.c_str() );
	if (infile.is_open())
	{
		string line;
		while ( getline (infile, line) )
		{
			istringstream iss(line);
			if ( linecount++ < 1 )						// skip first row
				continue;
			else
				iss >> dummy >> tau >> x >> y >> pin
					>> Tin >> muBin >> muSin >> muQin
					>> ein >> Bin >> Sin >> Qin >> sin;	// discard rest of line after this

			xvec.push_back( x );
			yvec.push_back( y );
			pvec.push_back( pin );
			Tvec.push_back( Tin );
			muBvec.push_back( muBin );
			muSvec.push_back( muSin );
			muQvec.push_back( muQin );
			evec.push_back( ein );
			Bvec.push_back( Bin );
			Svec.push_back( Sin );
			Qvec.push_back( Qin );
			svec.push_back( sin );
			nSPH++;
		}
	}

	infile.close();

	// now loop through and compute physical quantities
	const double dx = 0.1, dy = 0.1;
	const double xmin = -15.0, ymin = -15.0;
	const double xmax = -xmin, ymax = -ymin;

	ofstream outfile( outfilename.c_str() );

	for (double x_local = xmin; x_local < xmax + 1e-10; x_local += dx )
	for (double y_local = ymin; y_local < ymax + 1e-10; y_local += dy )
	{
		double normalization 				= 1e-100;	// protects from dividing by zero below
		double temperature 					= 0.0;
		double baryon_chemical_potential 	= 0.0;
		double strange_chemical_potential 	= 0.0;
		double electric_chemical_potential 	= 0.0;
		double energy_density 				= 0.0;
		double baryon_density 				= 0.0;
		double strange_density 				= 0.0;
		double electric_density 			= 0.0;
		double entropy_density 				= 0.0;

		// loop over SPH particles
		for (int iSPH = 0; iSPH < nSPH; iSPH++)
		{
			double kern 				 = kernel(x_local - xvec[iSPH], y_local - yvec[iSPH]);
			normalization 				+= kern;
			energy_density 				+= kern * evec[iSPH];
			baryon_density 				+= kern * Bvec[iSPH];
			strange_density 			+= kern * Svec[iSPH];
			electric_density 			+= kern * Qvec[iSPH];
			temperature 				+= kern * Tvec[iSPH];
			baryon_chemical_potential 	+= kern * muBvec[iSPH];
			strange_chemical_potential 	+= kern * muSvec[iSPH];
			electric_chemical_potential += kern * muQvec[iSPH];
			entropy_density 			+= kern * svec[iSPH];
		}

		// normalize results
		energy_density 				/= normalization;
		baryon_density 				/= normalization;
		strange_density 			/= normalization;
		electric_density 			/= normalization;
		temperature 				/= normalization;
		baryon_chemical_potential 	/= normalization;
		strange_chemical_potential 	/= normalization;
		electric_chemical_potential /= normalization;
		entropy_density 			/= normalization;

		outfile << setw(12) << setprecision(8) << scientific
			<< tau << "   " << x_local << "   " << y_local << "   "
			<< temperature << "   " << baryon_chemical_potential << "   "
			<< strange_chemical_potential << "   " << electric_chemical_potential << "   "
			<< energy_density << "   " << baryon_density << "   " << strange_density << "   "
			<< electric_density << "   " << entropy_density << endl;

	}

	outfile.close();

	return 0;
}
