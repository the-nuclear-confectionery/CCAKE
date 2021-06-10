#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string>

using namespace std;

double get_min_max( vector<double> & v, const size_t nTpts, const size_t nmuBpts,
										const size_t nmuSpts, const size_t nmuQpts )
{
	vector<double> maxima;
	size_t count = 0;
	for ( size_t iT = 0; iT < nTpts; iT++ )
	for ( size_t imuB = 0; imuB < nmuBpts; imuB++ )
	for ( size_t imuS = 0; imuS < nmuSpts; imuS++ )
	{
		maxima.push_back(
			*std::max_element( v.begin() + count*nmuQpts,
							   v.begin() + (count+1)*nmuQpts,
								[](const double & a, const double & b)
								{ return abs(a) < abs(b); }
							)
						);
		count++;
	}

	return ( *std::min_element( maxima.begin(), maxima.end() ) );

}
