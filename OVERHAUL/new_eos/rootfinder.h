#ifndef ROOTFINDER_H
#define ROOTFINDER_H

#include <fstream>
#include <iostream>
#include <vector>

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <lib.h>
#include "/projects/jnorhos/BSQ/OVERHAUL/constants.h"

using namespace std;
using namespace constants;

class Rootfinder
{
private:

  vector<double> tbqsPosition;

  void tbqs( vector<double> & tbqsIn );
  void tbqs(double setT, double setmuB, double setmuQ, double setmuS);

  int rootfinder_fsbqs(const gsl_vector *x, void *params, gsl_vector *f);
  int rootfinder_febqs(const gsl_vector *x, void *params, gsl_vector *f);

  bool rootfinder4D(double e_or_s_Given, int e_or_s_mode,
						double rhoBGiven, double rhoSGiven, double rhoQGiven,
						double error, size_t steps);

public:

  bool update_s( double sin, double Bin, double Sin, double Qin,
                 vector<double> & updated_tbqs );

  double s_out( double ein, double Bin, double Sin, double Qin,
                vector<double> & updated_tbqs );

  

};


#endif