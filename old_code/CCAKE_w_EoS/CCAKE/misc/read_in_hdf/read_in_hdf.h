#ifndef READ_IN_HDF_H
#define READ_IN_HDF_H

#include <iostream>
#include <string>
#include <vector>

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

using namespace std;

inline void read_in_hdf(vector<vector<double> > & v, string filename)
{
    const H5std_string	FILE_NAME(filename.c_str());
	const H5std_string	DATASET_NAME("EOS");

    try
    {
        Exception::dontPrint();

		hid_t fileID, dset;
		fileID = H5Fopen (filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		dset = H5Dopen (fileID, "EOS", H5P_DEFAULT);

        H5File file(FILE_NAME, H5F_ACC_RDONLY);
        DataSet dataset = file.openDataSet(DATASET_NAME);

		// find dimensions
		hid_t dspace = H5Dget_space(dset);
		const int ndims = H5Sget_simple_extent_ndims(dspace);

		hsize_t dims[ndims];
		H5Sget_simple_extent_dims(dspace, dims, NULL);

		//const long long		DIM0 = 274591955;
		//const long long		DIM1 = 11;
		const long long DIM0 = dims[0];
		const long long DIM1 = dims[1];

std::cout << "DIM0 = " << DIM0 << std::endl;
std::cout << "DIM1 = " << DIM1 << std::endl;

		double data[DIM0][DIM1];

        dataset.read(data, PredType::NATIVE_DOUBLE);

		// transfer to vector
		v.clear();
		v = vector<vector<double> >( DIM0, vector<double>( DIM1, 0.0 ) );

		for (long long ix = 0; ix < DIM0; ix++)
		for (long long iy = 0; iy < DIM1; iy++)
			v[ix][iy] = data[ix][iy];

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

    return;  // successfully terminated
}

#endif
