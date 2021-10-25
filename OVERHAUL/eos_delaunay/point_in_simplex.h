#ifndef POINT_IN_SIMPLEX_H
#define POINT_IN_SIMPLEX_H

#include <cmath>
#include <iostream>
#include <vector>

#include <gsl/gsl_linalg.h>

using namespace std;

void get_barycentric_coordinates(double T[], double d[], double lambda[], int dim);
bool point_is_in_simplex( const vector<vector<double> > & v, const vector<double> & p, 
							vector<double> & lambda, bool verbose = true );

#endif
