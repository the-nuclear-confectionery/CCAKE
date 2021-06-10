#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string>

using namespace std;

double get_min_max( vector<double> & v, const size_t nTpts, const size_t nmuBpts,
										const size_t nmuQpts, const size_t nmuSpts )
{
	vector<double> maxima;
	size_t count = 0;
	for ( size_t iT = 0; iT < nTpts; iT++ )
	for ( size_t imuB = 0; imuB < nmuBpts; imuB++ )
	for ( size_t imuQ = 0; imuQ < nmuQpts; imuQ++ )
	{
		maxima.push_back( abs(
			*std::max_element( v.begin() + count*nmuSpts,
							   v.begin() + (count+1)*nmuSpts,
								[](const double & a, const double & b)
								{ return abs(a) < abs(b); }
							) )
						);
		count++;
	}

	return ( *std::min_element( maxima.begin(), maxima.end() ) );

}
