#include <iostream>
#include <string>
#include <vector>

#include "read_in_hdf.h"

using namespace std;

int main (int argc, char ** argv)
{
	if ( argc != 2 ) exit(8);
    
	string filename = argv[1];

    vector<vector<double> > data;

	cout << "Starting..." << endl;
	read_in_hdf(data, filename);

	cout << "Finished!" << endl;

    return 0;
}
