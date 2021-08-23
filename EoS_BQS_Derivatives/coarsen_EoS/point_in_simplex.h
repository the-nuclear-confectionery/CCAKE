#ifndef POINT_IN_SIMPLEX_H
#define POINT_IN_SIMPLEX_H

#include <cmath>
#include <iostream>
#include <vector>

#include <gsl/gsl_linalg.h>

using namespace std;

void get_barycentric_coordinates(double T[], double d[], double lambda[], int dim);

void point_is_in_simplex( const vector<vector<double> > & v, const vector<double> & p )
{
	//const int dim = 3;
	const int dim = p.size();
	const int nVertices = dim + 1;	// by definition for a simplex

//	const double one_by_sqrt2 = 1.0/sqrt(2.0);
//	vector<double> v1 {2, 0, one_by_sqrt2};
//	vector<double> v2 {0, 1, -one_by_sqrt2};
//	vector<double> v3 {-1, 0, one_by_sqrt2};
//	vector<double> v4 {0, -1, -one_by_sqrt2};

//	vector<double> p {0,0,0};	// this is the point we want to check
//
//	vector<vector<double> > v;
//	v.push_back( v1 );
//	v.push_back( v2 );
//	v.push_back( v3 );
//	v.push_back( v4 );

	// construct T matrix
	double T[dim*dim];
	for (int i = 0; i < dim; i++)	// x, y, z, ...
	for (int j = 0; j < dim; j++)	// v1, v2, v3, ...
		T[i*dim+j] = v[j][i] - v[dim][i];

	// construct displacement vector
	double d[dim];
	for (int i = 0; i < dim; i++)	// x, y, z, ...
		d[i] = p[i] - v[dim][i];

	// construct lambda vector
	double lambda[dim];
	get_barycentric_coordinates(T, d, lambda, dim);

	bool point_in_simplex = true;

	// print results
	double lambda_sum = 0.0;
	for (int i = 0; i < dim; i++)
	{
		point_in_simplex = point_in_simplex && lambda[i] <= 1.0 && lambda[i] >= 0.0;
		lambda_sum += lambda[i];
		cout << "lambda[" << i << "] = " << lambda[i] << endl;
	}
	cout << "lambda[" << dim+1 << "] = " << 1.0 - lambda_sum << endl;
	point_in_simplex = point_in_simplex && lambda_sum <= 1.0 && lambda_sum >= 0.0;

	if (point_in_simplex)
		cout << "Point is in simplex!" << endl;
	else
		cout << "Point is not in simplex!" << endl;

	return 0;
}

void get_barycentric_coordinates(double T[], double d[], double lambda[], int dim)
{
	gsl_matrix_view m = gsl_matrix_view_array (T, dim, dim);
	gsl_vector_view b = gsl_vector_view_array (d, dim);
	gsl_vector *x = gsl_vector_alloc (dim);
	
	int s;
	gsl_permutation * p = gsl_permutation_alloc (dim);
	
	gsl_linalg_LU_decomp (&m.matrix, p, &s);
	gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

	for (int i = 0; i < dim; i++) lambda[i] = gsl_vector_get(x, i);
	
	gsl_permutation_free (p);
	gsl_vector_free (x);

	return;
}

#endif
