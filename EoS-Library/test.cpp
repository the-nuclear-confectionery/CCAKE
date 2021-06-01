#include <iostream>

#include "eos.h"

using namespace std;

int main( int argc, char ** argv )
{
	cout << "Hello, world!" << endl;

	bool using_HDF = true;
	eos EOS("inputfiles/quantityFile.h5", "inputfiles/derivFile.h5", 1, using_HDF);

	return 0;
}
