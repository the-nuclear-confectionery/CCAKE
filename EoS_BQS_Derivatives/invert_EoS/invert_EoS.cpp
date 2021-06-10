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
	if (argc < 6) exit(-1);

	// read path to input file from command line
	string path_to_file = string(argv[1]);
	const size_t nTpts = stoi(argv[2]);
	const size_t nmuBpts = stoi(argv[3]);
	const size_t nmuQpts = stoi(argv[4]);
	const size_t nmuSpts = stoi(argv[5]);

	vector<double> Tvec, muBvec, muSvec, muQvec;
	vector<double> evec, bvec, svec, qvec;

	// then read in file itself
	ifstream infile(path_to_file.c_str());
	if (infile.is_open())
	{
		size_t count = 0;
		string line;
		double dummy, Tin, muBin, muSin, muQin, bin, sin, qin, ein; 
		while ( getline (infile, line)
				and count < nTpts /*second condition just during debugging*/ )
		{
			istringstream iss(line);
			iss >> Tin >> muBin >> muQin >> muSin >> dummy >> dummy
				>> bin >> sin >> qin >> ein >> dummy;

			Tvec.push_back( Tin );
			muBvec.push_back( muBin );
			muSvec.push_back( muSin );
			muQvec.push_back( muQin );
			bvec.push_back( bin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) );				// 1/fm^3
			svec.push_back( sin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) );				// 1/fm^3
			qvec.push_back( qin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) );				// 1/fm^3
			evec.push_back( 0.001*ein*Tin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) );	// GeV/fm^3

			if (++count % 1000000 == 0) cout << "Read in " << count << " lines." << endl;
		}
	}

	infile.close();

cout << Tvec.size() << endl;

	// loop over each T value to get 
	for ( size_t iT = 0; iT < nTpts; iT++ )
	{
		// during debugging only
		if ( iT > 0 ) break;

		//============================
		// get min-max's in each b,q,s direction
	
		//============================
		double b_min_max = get_b_min_max(bvec, iT, nTpts, nmuBpts, nmuQpts, nmuSpts);
	
		cout << "b_min_max = " << Tvec[indexer( iT, 0, 0, 0, nTpts, nmuBpts, nmuQpts, nmuSpts )]
				<< "   " << b_min_max << endl;
		
		for (size_t i = 0; i < bvec.size(); i++)
			if ( abs( b_min_max - abs(bvec[i]) ) < 1e-10 )
				cout << "b range:" << i << "   " << Tvec[i] << "   " << muBvec[i] << "   "
					<< muSvec[i] << "   " << muQvec[i] << "   "
					<< evec[i] << "   " << bvec[i] << "   "
					<< svec[i] << "   " << qvec[i] << endl;
	
		//============================
		double s_min_max = get_s_min_max(svec, iT, nTpts, nmuBpts, nmuQpts, nmuSpts);
	
		cout << "s_min_max = " << Tvec[indexer( iT, 0, 0, 0, nTpts, nmuBpts, nmuQpts, nmuSpts )]
				<< "   " << s_min_max << endl;
		
		for (size_t i = 0; i < svec.size(); i++)
			if ( abs( s_min_max - abs(svec[i]) ) < 1e-10 )
				cout << "s range:" << i << "   " << Tvec[i] << "   " << muBvec[i] << "   "
					<< muSvec[i] << "   " << muQvec[i] << "   "
					<< evec[i] << "   " << bvec[i] << "   "
					<< svec[i] << "   " << qvec[i] << endl;
	
		//============================
		double q_min_max = get_q_min_max(qvec, iT, nTpts, nmuBpts, nmuQpts, nmuSpts);
	
		cout << "q_min_max = " << Tvec[indexer( iT, 0, 0, 0, nTpts, nmuBpts, nmuQpts, nmuSpts )]
				<< "   " << q_min_max << endl;
		
		for (size_t i = 0; i < qvec.size(); i++)
			if ( abs( q_min_max - abs(qvec[i]) ) < 1e-10 )
				cout << "q range:" << i << "   " << Tvec[i] << "   " << muBvec[i] << "   "
					<< muSvec[i] << "   " << muQvec[i] << "   "
					<< evec[i] << "   " << bvec[i] << "   "
					<< svec[i] << "   " << qvec[i] << endl;

		cout << endl << endl;
	}



	return 0;
}

