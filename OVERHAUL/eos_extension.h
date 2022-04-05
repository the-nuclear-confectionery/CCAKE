#ifndef EOS_EXTENSION_H
#define EOS_EXTENSION_H

#include <algorithm>
#include <cmath>
#include <iostream>

#include "constants.h"

namespace eos_extension
{
  using namespace constants;
  using std::cout;
  using std::endl;

	const double hc = hbarc_MeVfm;

  extern double b0000, b0002, b0020, b0200, b2000;
  extern double b0022, b0202, b2002, b0220, b2020, b2200;
  extern double b0004, b0040, b0400, b4000;

  //////////////////////////////////////////////////////////////////////////////
  // function definitions
  void set_coeffs( const double point[], double thermodynamics[] );

  double p(double T, double muB, double muQ, double muS);
	double s(double T, double muB, double muQ, double muS);
	double B(double T, double muB, double muQ, double muS);
	double S(double T, double muB, double muQ, double muS);
	double Q(double T, double muB, double muQ, double muS);
	double P2T2(double T, double muB, double muQ, double muS);
	double P2B2(double T, double muB, double muQ, double muS);
	double P2S2(double T, double muB, double muQ, double muS);
	double P2Q2(double T, double muB, double muQ, double muS);
	double P2BQ(double T, double muB, double muQ, double muS);
	double P2BS(double T, double muB, double muQ, double muS);
	double P2QS(double T, double muB, double muQ, double muS);
	double P2TB(double T, double muB, double muQ, double muS);
	double P2TS(double T, double muB, double muQ, double muS);
	double P2TQ(double T, double muB, double muQ, double muS);

  void get_full_thermo( const double point[], double results[] );
  void get_eBSQ( const double point[], double results[] );
  void get_sBSQ( const double point[], double results[] );

  void get_nonconformal_extension( const double point[], double thermodynamics[] );
  void project_to_boundary( double point[], const double minima[], const double maxima[] );
}

#endif