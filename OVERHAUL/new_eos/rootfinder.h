#ifndef ROOTFINDER_H
#define ROOTFINDER_H

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

class Rootfinder
{
private:
  //struct to pass the target (E, rhoB, rhoQ, rhoS) into the rootfinder function
  struct rootfinder_parameters
  {
      double eorEntGiven;          //these are the desired s and BSQ
      double rhoBGiven;
      double rhoQGiven;
      double rhoSGiven;
      rootfinder_parameters();
      rootfinder_parameters(
      double seteorEntGiven, double setRhoBGiven,
      double setRhoQGiven, double setRhoSGiven);
  public:
      void set( double setEorEntGiven, double setRhoBGiven,
          double setRhoQGiven, double setRhoSGiven);
  };

  //Default constructor to make the compiler happy. Should never be called
  rootfinder_parameters::rootfinder_parameters() {}
  //constructor which initializes all struct variables
  rootfinder_parameters::rootfinder_parameters(
    double setEorEntGiven, double setRhoBGiven,
    double setRhoQGiven, double setRhoSGiven
    )
  {
      eorEntGiven = setEorEntGiven;
      rhoBGiven = setRhoBGiven;
      rhoQGiven = setRhoQGiven;
      rhoSGiven = setRhoSGiven;
  }

  void rootfinder_parameters::set(
    double setEorEntGiven, double setRhoBGiven,
    double setRhoQGiven, double setRhoSGiven)
  {
      eorEntGiven = setEorEntGiven;
      rhoBGiven = setRhoBGiven;
      rhoQGiven = setRhoQGiven;
      rhoSGiven = setRhoSGiven;

  }

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