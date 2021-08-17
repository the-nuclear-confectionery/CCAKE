#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string>

using namespace std;

constexpr int mode = 1;
constexpr double EPS = 1e-25;
const double hbarc = 197.327;

bool are_nearby(const vector<double> & d1, const vector<double> & d2, double r0)
{
	if (mode == 1)
		return (  (d1[0]-d2[0])*(d1[0]-d2[0])
			+ (d1[1]-d2[1])*(d1[1]-d2[1])
			+ (d1[2]-d2[2])*(d1[2]-d2[2])
			+ (d1[3]-d2[3])*(d1[3]-d2[3]) < r0*r0 );
	else if (mode == 2)
		return (   (d1[0]-d2[0])*(d1[0]-d2[0]) < r0*r0
			&& (d1[1]-d2[1])*(d1[1]-d2[1]) < r0*r0
			&& (d1[2]-d2[2])*(d1[2]-d2[2]) < r0*r0
			&& (d1[3]-d2[3])*(d1[3]-d2[3]) < r0*r0 );
	else
                return (   (d1[0]-d2[0])*(d1[0]-d2[0]) < r0*r0
                        || (d1[1]-d2[1])*(d1[1]-d2[1]) < r0*r0
                        || (d1[2]-d2[2])*(d1[2]-d2[2]) < r0*r0
                        || (d1[3]-d2[3])*(d1[3]-d2[3]) < r0*r0 );

}


bool density_comp( const vector<double> & d1, const vector<double> & d2 )
{
	if (d1[0]==d2[0])
	{
		if (d1[1]==d2[1])
		{
			if (d1[2]==d2[2])
			{
				if (d1[3]==d2[3]) return true;
				else return (d1[3] < d2[3]);
			}
			else return (d1[2] < d2[2]);
		}
		else return (d1[1] < d2[1]);
	}
	else return (d1[0] < d2[0]);
}

int main(int argc, char *argv[])
{
	// check input first
	if (argc < 2) exit(-1);

	// read path to input file from command line
	string path_to_file = string(argv[1]);

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
			count++;
			//if (count % 100 != 0) continue;

			istringstream iss(line);
			iss >> Tin >> muBin >> muQin >> muSin >> dummy >> dummy
				>> bin >> sin >> qin >> ein >> dummy;

			Tvec.push_back( Tin );
			muBvec.push_back( muBin );
			muSvec.push_back( muSin );
			muQvec.push_back( muQin );
			bvec.push_back( bin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) );		// 1/fm^3
			svec.push_back( sin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) );		// 1/fm^3
			qvec.push_back( qin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) );		// 1/fm^3
			evec.push_back( ein*Tin*Tin*Tin*Tin/(hbarc*hbarc*hbarc) );	// MeV/fm^3

			if (count % 1000000 == 0) cout << "Read in " << count << " lines." << endl;
		}
	}

	infile.close();

	cout << "Initialize densities..." << endl;
	vector<double> logevec(evec.size()), rhoBtvec(bvec.size()),
			rhoStvec(svec.size()), rhoQtvec(qvec.size());
	vector<vector<double> > densities;
	for ( size_t ivec = 0; ivec < evec.size(); ivec++ )
	{
		logevec[ivec] = log(evec[ivec]);
		rhoBtvec[ivec] = bvec[ivec] / pow(evec[ivec]/hbarc, 0.75);
                rhoStvec[ivec] = svec[ivec] / pow(evec[ivec]/hbarc, 0.75);
                rhoQtvec[ivec] = qvec[ivec] / pow(evec[ivec]/hbarc, 0.75);
		vector<double> density(5); // 5th element is cell index
		density[0] = logevec[ivec];
		density[1] = rhoBtvec[ivec];
                density[2] = rhoStvec[ivec];
                density[3] = rhoQtvec[ivec];
		density[4] = ivec;
		densities.push_back( density );
	}

	double minloge = *min_element(logevec.begin(), logevec.end());
	double maxloge = *max_element(logevec.begin(), logevec.end());
	double minrBt = *min_element(rhoBtvec.begin(), rhoBtvec.end());
	double maxrBt = *max_element(rhoBtvec.begin(), rhoBtvec.end());
	double minrSt = *min_element(rhoStvec.begin(), rhoStvec.end());
	double maxrSt = *max_element(rhoStvec.begin(), rhoStvec.end());
	double minrQt = *min_element(rhoQtvec.begin(), rhoQtvec.end());
	double maxrQt = *max_element(rhoQtvec.begin(), rhoQtvec.end());

	// normalize to unit hypercube
	cout << "Normalizing..." << endl;
	for ( auto & density : densities )
	{
		density[0] = (density[0] - minloge) / ( maxloge - minloge );
                density[1] = (density[1] - minrBt) / ( maxrBt - minrBt );
                density[2] = (density[2] - minrSt) / ( maxrSt - minrSt );
                density[3] = (density[3] - minrQt) / ( maxrQt - minrQt );
	}

	cout << "Sorting..." << endl;
	std::sort(densities.begin(), densities.end(), density_comp);


	vector<vector<double> > coarsened_densities;

	cout << "Size before: " << densities.size() << endl;

	size_t n_deleted_elements = 0;

	//constexpr double r0 = 0.01;	// units normalized to unit hypercube
	const double r0 = ( argc > 2 ) ? stod(argv[2]) : 0.025;
	cout << "r0 = " << r0 << endl;
	/*while ( densities.size() > 0 )
	{
		// push a density to coarsened_densities and remove it from densities
		vector<double> current_density = densities.back();
		coarsened_densities.push_back( current_density );
		densities.pop_back();

		// loop through the rest of densities and remove neighbors inside norm<r
		size_t current_size = densities.size();
		cout << "densities.size() = " << current_size << endl;
		if (current_size == 0) break; // nothing left to do
		for (size_t i = current_size-1; i >= 0; i--)
		{
			if ( are_nearby(current_density, densities[i], r0) )
			{
				densities[i] = densities.back(); // overwrite with last element
				densities.pop_back();		 // remove last element
				n_deleted_elements++;
			}
			if (i==0) break; // need this since otherwise i loops to large +ve value
		}
	}*/

	int count = 0;
	for ( const auto & current_density : densities  )
	{
		if ( count%1000 == 0) cout << "count = " << count++ << endl;
		bool found_point_too_close = false;
		for ( const auto & cd : coarsened_densities )
		{
			if ( are_nearby(current_density, cd, r0) )	// if this density too close to
			{											// some coarsened density...
				found_point_too_close = true;
				break;
			}
		}

		if ( !found_point_too_close )
			coarsened_densities.push_back( current_density );
	}

	cout << "Size after: " << densities.size() << endl;
	cout << "coarsened_densities.size() = " << coarsened_densities.size() << endl;
	cout << "Deleted " << n_deleted_elements << " elements." << endl << endl;

	// sort coarsened density vector by (last) index
	cout << "Re-sorting..." << endl << endl;
	std::sort(coarsened_densities.begin(), coarsened_densities.end(),
			[](const vector<double> & a, const vector<double> & b)
			{ return a[4] < b[4]; } );

	cout << "Results:" << endl;
	for ( const auto & cd : coarsened_densities )
	{
		size_t idx = (size_t)cd.back();
		cout << idx << "   "
			<< Tvec[idx] << "   " << muBvec[idx] << "   "
                        << muSvec[idx] << "   " << muQvec[idx] << "   "
			<< cd[0] << "   " << cd[1] << "   " << cd[2] << "   " << cd[3] << endl;
	}

	return 0;
}

