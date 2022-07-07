#include <bitset>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

vector<int> to_binary_digits( int i )
{
	//const int n = static_cast<int>( ceil( log2(i) ) );
	bitset<8> bs(i);
	vector<int> digits;
	for (int iDigit = bs.size()-1; iDigit >= 0; iDigit--)
		digits.push_back( bs[iDigit] );
	return digits;
}

int main( int argc, char ** argv )
{
	cout << argv[1] << ": ";
	for ( auto & digit : to_binary_digits( stoi( argv[1] ) ) )
		cout << " " << digit;
	cout << endl;
	return 0;
}
