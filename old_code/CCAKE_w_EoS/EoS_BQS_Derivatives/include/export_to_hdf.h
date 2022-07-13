#ifndef EXPORT_TO_HDF_H
#define EXPORT_TO_HDF_H

#include "hdf5.h"

void export_to_HDF(double ** array, char filename[], long long gridLength, long long gridWidth)
{
	hid_t   file_id, dataset_id, dataspace_id;
	hsize_t dims[2];
	herr_t  status;

	double arrayToStore[gridLength][gridWidth];
	for (long long i = 0; i < gridLength; i+=1)
	for (long long j = 0; j < gridWidth;  j+=1)
		arrayToStore[i][j] = array[i][j];

	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	dims[0]      = gridLength;
	dims[1]      = gridWidth;
	dataspace_id = H5Screate_simple(2, dims, NULL);

	dataset_id = H5Dcreate2(file_id, "/EOS", H5T_NATIVE_DOUBLE, dataspace_id, 
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arrayToStore);

	status = H5Dclose(dataset_id);

	status = H5Sclose(dataspace_id);

  status = H5Fclose(file_id);
}


// use this to export full thermodynamics with grid dimensions along each axis
void export_thermo_to_HDF( double ** array, char filename[],
                           long long gridLength, long long gridWidth,
                           int grid_dimensions[4] )
{
	hid_t   file_id, dataset_id, dataspace_id;
	hsize_t grid_dims[1], dims[2];
	herr_t  status;

	double arrayToStore[gridLength][gridWidth];
	for (long long i = 0; i < gridLength; i+=1)
	for (long long j = 0; j < gridWidth;  j+=1)
		arrayToStore[i][j] = array[i][j];

	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  //------------------------------------
  // set all thermodynamics quantities
	dims[0]      = gridLength;
	dims[1]      = gridWidth;
	dataspace_id = H5Screate_simple(2, dims, NULL);

	dataset_id = H5Dcreate2(file_id, "/data", H5T_NATIVE_DOUBLE, dataspace_id, 
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arrayToStore);

	status = H5Dclose(dataset_id);

	status = H5Sclose(dataspace_id);
  //------------------------------------


  //------------------------------------
  // set grid dimensions
  grid_dims[0] = 4; // 4D phase diagram
	dataspace_id = H5Screate_simple(1, grid_dims, NULL);

	dataset_id = H5Dcreate2(file_id, "/dimensions", H5T_NATIVE_INT, dataspace_id, 
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, grid_dimensions);

	status = H5Dclose(dataset_id);

	status = H5Sclose(dataspace_id);
  //------------------------------------

  status = H5Fclose(file_id);
}






// develop this version later in order to store quantities with names, units, etc.
/*void export_thermo_to_HDF(double ** array, char filename[], long long gridLength, long long gridWidth)
{
	hid_t   file_id, dataset_id, dataspace_id;
	hsize_t dims[1];
	herr_t  status;

  // set filename to store data
	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // start with just the temperature
	double arrayToStore[gridLength];
	for (long long i = 0; i < gridLength; i+=1)
		arrayToStore[i][0] = array[i][0];


	dims[0]      = gridLength;
	dataspace_id = H5Screate_simple(1, dims, NULL);

	dataset_id = H5Dcreate2(file_id, "/EOS/T", H5T_NATIVE_DOUBLE, dataspace_id, 
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arrayToStore);



	status = H5Dclose(dataset_id);

	status = H5Sclose(dataspace_id);

  status = H5Fclose(file_id);
}*/

#endif
