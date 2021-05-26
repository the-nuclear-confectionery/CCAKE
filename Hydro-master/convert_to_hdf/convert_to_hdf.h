#ifndef CONVERT_TO_HDF_H
#define CONVERT_TO_HDF_H

#include <fstream>
#include <sstream>
#include <vector>

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

using namespace std;

void read_in_data(vector<vector<double> > & data, string filename, int nHeaderLines = 0)
{
	data.clear();

	ifstream infile(filename.c_str());

	if (infile.is_open())
	{
		string line;
		int count = 0;
		while ( getline (infile, line) )
		{
			if ( count++ < nHeaderLines ) continue;

			istringstream iss(line);
			vector<double> datum;
			double tmp;
			while ( iss >> tmp ) datum.push_back( tmp );
			
			data.push_back( datum );
		}
	}

	std::cout << "Finished reading in dataset of size = " << data.size() << std::endl;

	infile.close();
	return;
}

void output_to_HDF( const vector<vector<double> > & v, string outfilename )
{
	std::cout << "Writing dataset out to " << outfilename << std::endl;

	const H5std_string	FILE_NAME(outfilename.c_str());
	const H5std_string	DATASET_NAME("EOS");
	const long long	 	NX = v.size();
	const long long	 	NY = (v[0]).size();
	const long long	 	RANK = 2;

	double data[NX][NY];
	for (long long ix = 0; ix < NX; ix++)
	for (long long iy = 0; iy < NY; iy++)
		data[ix][iy] = v[ix][iy];

    try
    {
		Exception::dontPrint();
	
		H5File file(FILE_NAME, H5F_ACC_TRUNC);
	
		hsize_t dims[2];
		dims[0] = NX;
		dims[1] = NY;
		DataSpace dataspace(RANK, dims);
	
		DataSet dataset = file.createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, dataspace);
			
		dataset.write(data, PredType::NATIVE_DOUBLE);
    }

    catch(FileIException error)
    {
		error.printError();
		return;
    }

    catch(DataSetIException error)
    {
		error.printError();
		return;
    }

    catch(DataSpaceIException error)
    {
		error.printError();
		return;
    }
}

#endif
