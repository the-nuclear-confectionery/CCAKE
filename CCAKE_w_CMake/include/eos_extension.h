#ifndef EOS_EXTENSION_H
#define EOS_EXTENSION_H

#include <algorithm>
#include <cmath>
#include <iostream>

#include "constants.h"
#include "eos_header.h"

namespace eos_extension
{
  using namespace constants;
  using std::cout;
  using std::endl;

	const double hc = hbarc_MeVfm;

  extern double b0000, b0001, b0010, b0100, b1000;
  extern double b0011, b0101, b1001, b0110, b1010, b1100;
  extern double b0002, b0020, b0200, b2000;

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

  void get_full_thermo( const double point[], double results[] ); // results length = 17
  void get_eBSQ( const double point[], double results[] );        // results length = 4
  void get_sBSQ( const double point[], double results[] );        // results length = 4

  // projects point on boundary to arbitrary point outside phase diagram boundary
  void get_nonconformal_extension( const double point[],
                                   const double point_projected[],
                                   double results[] );

  void project_to_boundary( double point[], const double minima[],
                                            const double maxima[] );
}

#endif