#ifndef TOY_THERMO_H
#define TOY_THERMO_H

#include "constants.h"

namespace toy_thermo
{
  using namespace constants;

	const double hc = hbarc_MeVfm;
	const double toy_T_scale = 150.0, toy_muB_scale = 101.0,
               toy_muS_scale = 103.0, toy_muQ_scale = 113.0;

	double p(double T, double muB, double muQ, double muS);
	double s(double T, double muB, double muQ, double muS);
	double B(double T, double muB, double muQ, double muS);
  double S(double T, double muB, double muQ, double muS);
	double Q(double T, double muB, double muQ, double muS);
	double P2B2(double T, double muB, double muQ, double muS);
	double P2Q2(double T, double muB, double muQ, double muS);
	double P2S2(double T, double muB, double muQ, double muS);
	double P2BQ(double T, double muB, double muQ, double muS);
	double P2BS(double T, double muB, double muQ, double muS);
	double P2QS(double T, double muB, double muQ, double muS);
	double P2TB(double T, double muB, double muQ, double muS);
	double P2TQ(double T, double muB, double muQ, double muS);
	double P2TS(double T, double muB, double muQ, double muS);
	double P2T2(double T, double muB, double muQ, double muS);

}


#endif