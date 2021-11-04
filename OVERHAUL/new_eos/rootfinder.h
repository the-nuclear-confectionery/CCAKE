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
public:
  Rootfinder(){}
  ~Rootfinder(){}

private:

  const int VERBOSE = 5;
  const int STEPS = 1000000;
  const double TOLERANCE = 1e-12;

  //Rootfinding method used **THIS CAN BE CHANGED DEPENDING ON EOS
  const gsl_multiroot_fsolver_type *TYPE = gsl_multiroot_fsolver_hybrids;

  double minT, maxT, minMuB, maxMuB, minMuS, maxMuS, minMuQ, maxMuQ;

  vector<double> tbqsPosition;

  void tbqs( vector<double> & tbqsIn );
  void tbqs(double setT, double setmuB, double setmuQ, double setmuS);

  int rootfinder_fsbqs(const gsl_vector *x, void *params, gsl_vector *f);
  int rootfinder_febqs(const gsl_vector *x, void *params, gsl_vector *f);

  bool rootfinder4D(double e_or_s_Given, int e_or_s_mode,
						double rhoBGiven, double rhoSGiven, double rhoQGiven,
						double error, size_t steps);

public:

  // formerly update_s and s_out, respectively.
  bool find_sBSQ_root( double sin, double Bin, double Sin, double Qin,
                       vector<double> & updated_tbqs );

  bool find_eBSQ_root( double ein, double Bin, double Sin, double Qin,
                       vector<double> & updated_tbqs );

  

};


#endif