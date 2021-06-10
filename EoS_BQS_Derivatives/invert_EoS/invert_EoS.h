#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <string>

using namespace std;

size_t indexer( size_t iT, size_t imuB, size_t imuQ, size_t imuS,
				const size_t nTpts, const size_t nmuBpts,
				const size_t nmuQpts, const size_t nmuSpts )
{
	return ( ( ( iT * nmuBpts + imuB ) * nmuQpts + imuQ ) * nmuSpts + imuS );
}


double get_e_min_max( vector<double> & v, const size_t nTpts, const size_t nmuBpts,
										const size_t nmuQpts, const size_t nmuSpts )
{
	vector<double> maxima;
	size_t count = 0;
	for ( size_t imuB = 0; imuB < nmuBpts; imuB++ )
	for ( size_t imuS = 0; imuS < nmuSpts; imuS++ )
	for ( size_t imuQ = 0; imuQ < nmuQpts; imuQ++ )
	{
		vector<double> slice;
		for ( size_t iT = 0; iT < nTpts; iT++ )
			slice.push_back( v[indexer( iT, imuB, imuQ, imuS, nTpts, nmuBpts, nmuQpts, nmuSpts )] );
		maxima.push_back( abs(
			*std::max_element( slice.begin(), slice.end(),
								[](const double & a, const double & b)
								{ return abs(a) < abs(b); }
							) )
						);
		count++;
	}

	return ( *std::min_element( maxima.begin(), maxima.end() ) );
}




double get_b_min_max( vector<double> & v, const size_t nTpts, const size_t nmuBpts,
										const size_t nmuQpts, const size_t nmuSpts )
{
	vector<double> maxima;
	size_t count = 0;
	for ( size_t iT = 0; iT < nTpts; iT++ )
	for ( size_t imuS = 0; imuS < nmuSpts; imuS++ )
	for ( size_t imuQ = 0; imuQ < nmuQpts; imuQ++ )
	{
		vector<double> slice;
		for ( size_t imuB = 0; imuB < nmuBpts; imuB++ )
			slice.push_back( v[indexer( iT, imuB, imuQ, imuS, nTpts, nmuBpts, nmuQpts, nmuSpts )] );
		maxima.push_back( abs(
			*std::max_element( slice.begin(), slice.end(),
								[](const double & a, const double & b)
								{ return abs(a) < abs(b); }
							) )
						);
		count++;
	}

	return ( *std::min_element( maxima.begin(), maxima.end() ) );
}


double get_q_min_max( vector<double> & v, const size_t nTpts, const size_t nmuBpts,
										const size_t nmuQpts, const size_t nmuSpts )
{
	vector<double> maxima;
	size_t count = 0;
	for ( size_t iT = 0; iT < nTpts; iT++ )
	for ( size_t imuB = 0; imuB < nmuBpts; imuB++ )
	for ( size_t imuS = 0; imuS < nmuSpts; imuS++ )
	{
		vector<double> slice;
		for ( size_t imuQ = 0; imuQ < nmuQpts; imuQ++ )
			slice.push_back( v[indexer( iT, imuB, imuQ, imuS, nTpts, nmuBpts, nmuQpts, nmuSpts )] );
		maxima.push_back( abs(
			*std::max_element( slice.begin(), slice.end(),
								[](const double & a, const double & b)
								{ return abs(a) < abs(b); }
							) )
						);
		count++;
	}

	return ( *std::min_element( maxima.begin(), maxima.end() ) );
}

double get_s_min_max( vector<double> & v, const size_t nTpts, const size_t nmuBpts,
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

