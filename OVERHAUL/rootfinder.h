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

#include <lib.h>
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

  const int VERBOSE = 10;
  const int STEPS = 1000;
  const double TOLERANCE = 1e-12;

  //Rootfinding method used **THIS CAN BE CHANGED DEPENDING ON EOS
  const gsl_multiroot_fsolver_type *TYPE = gsl_multiroot_fsolver_hybrids;

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

  /*void set_grid_ranges( double minT_in, double maxT_in, double minMuB_in, double maxMuB_in,
                        double minMuS_in, double maxMuS_in, double minMuQ_in, double maxMuQ_in )
  {
    minT   = minT_in;   maxT   = maxT_in;
    minMuB = minMuB_in; maxMuB = maxMuB_in;
    minMuS = minMuS_in; maxMuS = maxMuS_in;
    minMuQ = minMuQ_in; maxMuQ = maxMuQ_in;
  }*/

  // formerly update_s and s_out, respectively.
  bool find_sBSQ_root( double sin, double Bin, double Sin, double Qin,
                       std::function<void(double[], double[])> function_to_evaluate,
                       vector<double> & tbqs_minima,
                       vector<double> & tbqs_maxima,
                       vector<double> & updated_tbqs );

  bool find_eBSQ_root( double ein, double Bin, double Sin, double Qin,
                       std::function<void(double[], double[])> function_to_evaluate,
                       vector<double> & tbqs_minima,
                       vector<double> & tbqs_maxima,
                       vector<double> & updated_tbqs );

  

};


#endif