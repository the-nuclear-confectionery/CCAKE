#ifndef ROOTFINDER_H
#define ROOTFINDER_H

#include <fstream>
#include <functional>
#include <iostream>
#include <vector>

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "constants.h"

using namespace std;
using namespace constants;

int rootfinder_f( const gsl_vector *x, void *params, gsl_vector *f );


class Rootfinder
{
public:
  Rootfinder(){}
  ~Rootfinder(){}

private:

  const int VERBOSE = 0;
  const int STEPS = 1000;
  const double TOLERANCE = 1e-10;

  //Rootfinding method used **THIS CAN BE CHANGED DEPENDING ON EOS
//  const gsl_multiroot_fsolver_type *TYPE = gsl_multiroot_fsolver_hybrids;
//  const gsl_multiroot_fsolver_type *TYPE = gsl_multiroot_fsolver_hybrid;
//  const gsl_multiroot_fsolver_type *TYPE = gsl_multiroot_fsolver_dnewton;

  double minT, maxT, minMuB, maxMuB, minMuS, maxMuS, minMuQ, maxMuQ;

  vector<double> tbqsPosition;

  void tbqs( vector<double> & tbqsIn );
  void tbqs(double setT, double setmuB, double setmuQ, double setmuS);

  bool rootfinder4D(double e_or_s_Given, int e_or_s_mode,
						double rhoBGiven, double rhoSGiven, double rhoQGiven,
						double error, size_t steps,
            std::function<void(double[], double[])> function_to_evaluate,
            vector<double> & updated_tbqs );

public:
  
  bool find_root( const string & e_or_s, double ein_or_sin,
                  double Bin, double Sin, double Qin,
                  std::function<void(double[], double[])> function_to_evaluate,
                  vector<double> & tbqs_minima,
                  vector<double> & tbqs_maxima,
                  vector<double> & updated_tbqs );

};


#endif