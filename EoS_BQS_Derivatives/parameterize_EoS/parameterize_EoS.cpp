#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string>

#include "parameterize_EoS.h"

using namespace std;

const double hbarc = 197.327;

int main(int argc, char *argv[])
{
	// turn off errors for now...
	gsl_set_error_handler_off();

	// check input first
	if (argc >= 6)
	{
		// read path to input file from command line
		string path_to_file = string(argv[1]);
		const int nTpts = stoi(argv[2]);
		const int nmuBpts = stoi(argv[3]);
		const int nmuSpts = stoi(argv[4]);
		const int nmuQpts = stoi(argv[5]);

		vector<double> Tvec, muBvec, muSvec, muQvec;
		vector<double> evec, bvec, svec, qvec;
	
		// then read in file itself
		ifstream infile(path_to_file.c_str());
		if (infile.is_open())
		{
			int count = 0;
			string line;
			double dummy, Tin, muBin, muSin, muQin, bin, sin, qin, ein; 
			while ( getline (infile, line) )
			{
				istringstream iss(line);
				iss >> Tin >> muBin >> muQin >> muSin >> dummy >> dummy
					>> bin >> sin >> qin >> ein;
	
				Tvec.push_back( Tin );
				muBvec.push_back( muBin );
				muSvec.push_back( muSin );
				muQvec.push_back( muQin );
				bvec.push_back( bin );
				svec.push_back( sin );
				qvec.push_back( qin );
				evec.push_back( ein );

				if (++count % 100000 == 0) cout << "Read in " << count << " lines." << endl;
			}
		}
	
		infile.close();

		// just look at a slice for testing purposes
		vector<double> muBslice( muBvec.begin(), muBvec.begin() + nmuBpts*nmuSpts*nmuQpts );
		vector<double> muSslice( muSvec.begin(), muSvec.begin() + nmuBpts*nmuSpts*nmuQpts );
		vector<double> muQslice( muQvec.begin(), muQvec.begin() + nmuBpts*nmuSpts*nmuQpts );
		vector<double> eslice(   evec.begin(),   evec.begin()   + nmuBpts*nmuSpts*nmuQpts );
		vector<double> bslice(   bvec.begin(),   bvec.begin()   + nmuBpts*nmuSpts*nmuQpts );
		vector<double> sslice(   svec.begin(),   svec.begin()   + nmuBpts*nmuSpts*nmuQpts );
		vector<double> qslice(   qvec.begin(),   qvec.begin()   + nmuBpts*nmuSpts*nmuQpts );

		fit( muBslice, muSslice, muQslice, eslice, nmuBpts, nmuSpts, nmuQpts, true);
		fit( muBslice, muSslice, muQslice, bslice, nmuBpts, nmuSpts, nmuQpts, false);
		fit( muBslice, muSslice, muQslice, sslice, nmuBpts, nmuSpts, nmuQpts, false);
		fit( muBslice, muSslice, muQslice, qslice, nmuBpts, nmuSpts, nmuQpts, false);
	}

	return 0;
}





