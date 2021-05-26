#include <iostream>
#include <string>
#include <vector>

#include "convert_to_hdf.h"

using namespace std;

int main (int argc, char ** argv)
{
	if ( argc != 3 ) exit(8);

	string file_to_convert = argv[1];
	string converted_file  = argv[2];

	vector<vector<double> > data;

	// load file into dataset
	read_in_data(data, file_to_convert);

	// export dataset
	output_to_HDF(data, converted_file);

	return 0;
}
