#include <iostream>

#include "interpolatorND.h"

using namespace std;

int main()
{
	InterpolatorND<4> f;
	f.initialize( "data.dat" );

	vector<double> coords = {0.1,0.2,0.3,0.4};
	vector<double> results;
	f.evaluate( coords, results );

	cout << "results =" << endl;
	for ( auto & result : results ) cout << "   " << result << endl;

	InterpolatorND<2> f2;
	f2.initialize( "data2D.dat" );
	vector<double> coords2 = {0.1,0.2};
	vector<double> results2;

	f2.evaluate( coords2, results2 );
	cout << "results2 =" << endl;
	for ( auto & result : results2 ) cout << "   " << result << endl;

	return 0;
}
