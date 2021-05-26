#ifndef READ_IN_HDF_H
#define READ_IN_HDF_H

#include <iostream>
#include <string>

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

void read_in_hdf(vector<vector<double> > & v, string filename)
{
    const H5std_string	FILE_NAME(filename.c_str());
	const H5std_string	DATASET_NAME("EOS");
	const long long		DIM0 = 274591955;
	const long long		DIM1 = 11;

    long long i, j;
    double data[DIM0][DIM1];

	v.clear();
	v = vector<vector<double> >( DIM0, vector<double>( DIM1, 0.0 ) );

    try
    {
        Exception::dontPrint();

        H5File file(FILE_NAME, H5F_ACC_RDWR);
        DataSet dataset = file.openDataSet(DATASET_NAME);

        dataset.read(data, PredType::NATIVE_DOUBLE);

    }

	for (long long ix = 0; ix < NX; ix++)
	for (long long iy = 0; iy < NY; iy++)
		v[ix][iy] = data[ix][iy];


    catch(FileIException error)
    {
        error.printError();
        return -1;
    }

    // catch failure caused by the DataSet operations
    catch(DataSetIException error)
    {
        error.printError();
        return -1;
    }

    return 0;  // successfully terminated
}

#endif
