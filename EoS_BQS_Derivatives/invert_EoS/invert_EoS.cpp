#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string>

#include "invert_EoS.h"

using namespace std;

const double hbarc = 197.327;

int main(int argc, char *argv[])
{
	// check input first
	if (argc >= 6)
	{
		// read path to input file from command line
		string path_to_file = string(argv[1]);
		const int nTpts = stoi(argv[2]);
		const int nmuBpts = stoi(argv[3]);
		const int nmuQpts = stoi(argv[4]);
		const int nmuSpts = stoi(argv[5]);

		vector<double> Tvec, muBvec, muSvec, muQvec;
		vector<double> evec, bvec, svec, qvec;
	
		// then read in file itself
		ifstream infile(path_to_file.c_str());
		if (infile.is_open())
		{
			size_t count = 0;
			string line;
			double dummy, Tin, muBin, muSin, muQin, bin, sin, qin, ein; 
			while ( getline (infile, line) )
			{
				istringstream iss(line);
				iss >> Tin >> muBin >> muQin >> muSin >> dummy >> dummy
					>> bin >> sin >> qin >> ein >> dummy;
	
				Tvec.push_back( Tin );
				muBvec.push_back( muBin );
				muSvec.push_back( muSin );
				muQvec.push_back( muQin );
				bvec.push_back( bin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) );
				svec.push_back( sin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) );
				qvec.push_back( qin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) );
				evec.push_back( 0.001*ein*Tin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) );

				if (++count % 1000000 == 0) cout << "Read in " << count << " lines." << endl;
			}
		}
	
		infile.close();

		/*
		// just look at a slice for testing purposes
		vector<double> muBslice( muBvec.begin(), muBvec.begin() + nmuBpts*nmuSpts*nmuQpts );
		vector<double> muSslice( muSvec.begin(), muSvec.begin() + nmuBpts*nmuSpts*nmuQpts );
		vector<double> muQslice( muQvec.begin(), muQvec.begin() + nmuBpts*nmuSpts*nmuQpts );
		vector<double> eslice(   evec.begin(),   evec.begin()   + nmuBpts*nmuSpts*nmuQpts );
		vector<double> bslice(   bvec.begin(),   bvec.begin()   + nmuBpts*nmuSpts*nmuQpts );
		vector<double> sslice(   svec.begin(),   svec.begin()   + nmuBpts*nmuSpts*nmuQpts );
		vector<double> qslice(   qvec.begin(),   qvec.begin()   + nmuBpts*nmuSpts*nmuQpts );
		*/

		double s_min_max = get_min_max(svec, nTpts, nmuBpts, nmuQpts, nmuSpts);

		cout << "s_min_max = " << s_min_max << endl;
		
		for (size_t i = 0; i < svec.size(); i++)
			if ( abs( s_min_max - abs(svec[i]) ) < 1e-10 )
				cout << i << "   " << Tvec[i] << "   " << muBvec[i] << "   "
					<< muSvec[i] << "   " << muQvec[i] << "   "
					<< evec[i] << "   " << bvec[i] << "   "
					<< svec[i] << "   " << qvec[i] << endl;

	}

	return 0;
}





