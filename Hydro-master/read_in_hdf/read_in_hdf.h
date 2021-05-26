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

    catch(FileIException error)
    {
        error.printError();
        return;
    }

    // catch failure caused by the DataSet operations
    catch(DataSetIException error)
    {
        error.printError();
        return;
    }

	for (long long ix = 0; ix < DIM0; ix++)
	for (long long iy = 0; iy < DIM1; iy++)
		v[ix][iy] = data[ix][iy];

    return;  // successfully terminated
}

#endif
