#include "constants.h"
#include "toy_thermo.h"

namespace toy_thermo
{
  using namespace constants;

	double p(double T, double muB, double muQ, double muS)
	{
		double x = (T/toy_T_scale)*(T/toy_T_scale) + (muB/toy_muB_scale)*(muB/toy_muB_scale)
					 + (muS/toy_muS_scale)*(muS/toy_muS_scale) + (muQ/toy_muQ_scale)*(muQ/toy_muQ_scale);
		return x*x;
	}
	
	double s(double T, double muB, double muQ, double muS)
	{
		double x = (T/toy_T_scale)*(T/toy_T_scale) + (muB/toy_muB_scale)*(muB/toy_muB_scale)
					 + (muS/toy_muS_scale)*(muS/toy_muS_scale) + (muQ/toy_muQ_scale)*(muQ/toy_muQ_scale);
		return 4.0*hc*T*x/(toy_T_scale*toy_T_scale);
	}

	double B(double T, double muB, double muQ, double muS)
	{
		double x = (T/toy_T_scale)*(T/toy_T_scale) + (muB/toy_muB_scale)*(muB/toy_muB_scale)
					 + (muS/toy_muS_scale)*(muS/toy_muS_scale) + (muQ/toy_muQ_scale)*(muQ/toy_muQ_scale);
		return 4.0*hc*muB*x/(toy_muB_scale*toy_muB_scale);
	}

	double S(double T, double muB, double muQ, double muS)
	{
		double x = (T/toy_T_scale)*(T/toy_T_scale) + (muB/toy_muB_scale)*(muB/toy_muB_scale)
					 + (muS/toy_muS_scale)*(muS/toy_muS_scale) + (muQ/toy_muQ_scale)*(muQ/toy_muQ_scale);
		return 4.0*hc*muS*x/(toy_muS_scale*toy_muS_scale);
	}

	double Q(double T, double muB, double muQ, double muS)
	{
		double x = (T/toy_T_scale)*(T/toy_T_scale) + (muB/toy_muB_scale)*(muB/toy_muB_scale)
					 + (muS/toy_muS_scale)*(muS/toy_muS_scale) + (muQ/toy_muQ_scale)*(muQ/toy_muQ_scale);
		return 4.0*hc*muQ*x/(toy_muQ_scale*toy_muQ_scale);
	}

	double P2B2(double T, double muB, double muQ, double muS)
	{
		double x = (T/toy_T_scale)*(T/toy_T_scale) + (muB/toy_muB_scale)*(muB/toy_muB_scale)
					 + (muS/toy_muS_scale)*(muS/toy_muS_scale) + (muQ/toy_muQ_scale)*(muQ/toy_muQ_scale);
		return (8.0*muB*muB/(toy_muB_scale*toy_muB_scale*toy_muB_scale*toy_muB_scale)
				+ 4.0*x/(toy_muB_scale*toy_muB_scale))*hc*hc;
	}
	double P2Q2(double T, double muB, double muQ, double muS)
	{
		double x = (T/toy_T_scale)*(T/toy_T_scale) + (muB/toy_muB_scale)*(muB/toy_muB_scale)
					 + (muS/toy_muS_scale)*(muS/toy_muS_scale) + (muQ/toy_muQ_scale)*(muQ/toy_muQ_scale);
		return (8.0*muQ*muQ/(toy_muQ_scale*toy_muQ_scale*toy_muQ_scale*toy_muQ_scale)
				+ 4.0*x/(toy_muQ_scale*toy_muQ_scale))*hc*hc;
	}
	double P2S2(double T, double muB, double muQ, double muS)
	{
		double x = (T/toy_T_scale)*(T/toy_T_scale) + (muB/toy_muB_scale)*(muB/toy_muB_scale)
					 + (muS/toy_muS_scale)*(muS/toy_muS_scale) + (muQ/toy_muQ_scale)*(muQ/toy_muQ_scale);
		return (8.0*muS*muS/(toy_muS_scale*toy_muS_scale*toy_muS_scale*toy_muS_scale)
				+ 4.0*x/(toy_muS_scale*toy_muS_scale))*hc*hc;
	}
	
	double P2BQ(double T, double muB, double muQ, double muS)
			{ return (8.0*muB*muQ/(toy_muB_scale*toy_muB_scale*toy_muQ_scale*toy_muQ_scale))*hc*hc; }
	double P2BS(double T, double muB, double muQ, double muS)
			{ return (8.0*muB*muS/(toy_muB_scale*toy_muB_scale*toy_muS_scale*toy_muS_scale))*hc*hc; }
	double P2QS(double T, double muB, double muQ, double muS)
			{ return (8.0*muS*muQ/(toy_muS_scale*toy_muS_scale*toy_muQ_scale*toy_muQ_scale))*hc*hc; }
	
	double P2TB(double T, double muB, double muQ, double muS)
			{ return (8.0*T*muB/(toy_T_scale*toy_T_scale*toy_muB_scale*toy_muB_scale))*hc*hc; }
	double P2TQ(double T, double muB, double muQ, double muS)
			{ return (8.0*T*muQ/(toy_T_scale*toy_T_scale*toy_muQ_scale*toy_muQ_scale))*hc*hc; }
	double P2TS(double T, double muB, double muQ, double muS)
			{ return (8.0*T*muS/(toy_T_scale*toy_T_scale*toy_muS_scale*toy_muS_scale))*hc*hc; }
	double P2T2(double T, double muB, double muQ, double muS)
	{
		double x = (T/toy_T_scale)*(T/toy_T_scale) + (muB/toy_muB_scale)*(muB/toy_muB_scale)
					 + (muS/toy_muS_scale)*(muS/toy_muS_scale) + (muQ/toy_muQ_scale)*(muQ/toy_muQ_scale);
		return (8.0*T*T/(toy_T_scale*toy_T_scale*toy_T_scale*toy_T_scale)
				+ 4.0*x/(toy_T_scale*toy_T_scale))*hc*hc;
	}



}
