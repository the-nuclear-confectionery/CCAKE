#ifndef EOS_CONFORMAL_H
#define EOS_CONFORMAL_H

#include "constants.h"
#include "eos_base.h"
#include "eos_header.h"

class EoS_conformal: public EoS_base
{
private:
	const double hc = constants::hbarc_MeVfm;

  double c, T0, muB0, muS0, muQ0;

public:
  // default constructor/destructor
  EoS_conformal(){}
  ~EoS_conformal(){}

  { 0.0,     -TBQS_INFINITY, -TBQS_INFINITY, -TBQS_INFINITY };

  EoS_conformal( const double c_in, const double T0_in, const double muB0_in,
                 const double muS0_in, const double muQ0_in )
    { c = c_in; T0 = T0_in; muB0 = muB0_in; muS0 = muS0_in; muQ0 = muQ0_in; }


  double p(double T, double muB, double muQ, double muS)
	{
		double x = (T/T0)*(T/T0) + (muB/muB0)*(muB/muB0)
                + (muS/muS0)*(muS/muS0) + (muQ/muQ0)*(muQ/muQ0);
    double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
		return cp*x*x;
	}
	
	double s(double T, double muB, double muQ, double muS)
	{
		double x = (T/T0)*(T/T0) + (muB/muB0)*(muB/muB0)
                + (muS/muS0)*(muS/muS0) + (muQ/muQ0)*(muQ/muQ0);
    double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
		return 4.0*cp*hc*T*x/(T0*T0);
	}

	double B(double T, double muB, double muQ, double muS)
	{
		double x = (T/T0)*(T/T0) + (muB/muB0)*(muB/muB0)
                + (muS/muS0)*(muS/muS0) + (muQ/muQ0)*(muQ/muQ0);
    double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
		return 4.0*cp*hc*muB*x/(muB0*muB0);
	}

	double S(double T, double muB, double muQ, double muS)
	{
		double x = (T/T0)*(T/T0) + (muB/muB0)*(muB/muB0)
                + (muS/muS0)*(muS/muS0) + (muQ/muQ0)*(muQ/muQ0);
    double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
		return 4.0*cp*hc*muS*x/(muS0*muS0);
	}

	double Q(double T, double muB, double muQ, double muS)
	{
		double x = (T/T0)*(T/T0) + (muB/muB0)*(muB/muB0)
                + (muS/muS0)*(muS/muS0) + (muQ/muQ0)*(muQ/muQ0);
    double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
		return 4.0*cp*hc*muQ*x/(muQ0*muQ0);
	}

	double P2B2(double T, double muB, double muQ, double muS)
	{
		double x = (T/T0)*(T/T0) + (muB/muB0)*(muB/muB0)
                + (muS/muS0)*(muS/muS0) + (muQ/muQ0)*(muQ/muQ0);
    double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
		return cp*(8.0*muB*muB/(muB0*muB0*muB0*muB0) + 4.0*x/(muB0*muB0))*hc*hc;
	}
	double P2Q2(double T, double muB, double muQ, double muS)
	{
		double x = (T/T0)*(T/T0) + (muB/muB0)*(muB/muB0)
                + (muS/muS0)*(muS/muS0) + (muQ/muQ0)*(muQ/muQ0);
    double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
		return cp*(8.0*muQ*muQ/(muQ0*muQ0*muQ0*muQ0) + 4.0*x/(muQ0*muQ0))*hc*hc;
	}
	double P2S2(double T, double muB, double muQ, double muS)
	{
		double x = (T/T0)*(T/T0) + (muB/muB0)*(muB/muB0)
                + (muS/muS0)*(muS/muS0) + (muQ/muQ0)*(muQ/muQ0);
    double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
		return cp*(8.0*muS*muS/(muS0*muS0*muS0*muS0) + 4.0*x/(muS0*muS0))*hc*hc;
	}
	
	double P2BQ(double T, double muB, double muQ, double muS)
      { double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
        return (8.0*cp*muB*muQ/(muB0*muB0*muQ0*muQ0))*hc*hc; }
	double P2BS(double T, double muB, double muQ, double muS)
			{ double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
        return (8.0*cp*muB*muS/(muB0*muB0*muS0*muS0))*hc*hc; }
	double P2QS(double T, double muB, double muQ, double muS)
			{ double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
        return (8.0*cp*muS*muQ/(muS0*muS0*muQ0*muQ0))*hc*hc; }
	
	double P2TB(double T, double muB, double muQ, double muS)
			{ double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
        return (8.0*cp*T*muB/(T0*T0*muB0*muB0))*hc*hc; }
	double P2TQ(double T, double muB, double muQ, double muS)
			{ double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
        return (8.0*cp*T*muQ/(T0*T0*muQ0*muQ0))*hc*hc; }
	double P2TS(double T, double muB, double muQ, double muS)
			{ double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
        return (8.0*cp*T*muS/(T0*T0*muS0*muS0))*hc*hc; }
	double P2T2(double T, double muB, double muQ, double muS)
	{
    double x = (T/T0)*(T/T0) + (muB/muB0)*(muB/muB0)
                + (muS/muS0)*(muS/muS0) + (muQ/muQ0)*(muQ/muQ0);
    double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
		return cp*(8.0*T*T/(T0*T0*T0*T0) + 4.0*x/(T0*T0))*hc*hc;
	}

  void get_eBSQ( double point[], double results[] )
  {
    // point: (T, muB, muQ, muS)
    const double Tsol  = hc*point[0],
      muBsol           = hc*point[1],
      muQsol           = hc*point[2],
      muSsol           = hc*point[3];
    double POut        = p(Tsol, muBsol, muQsol, muSsol);
    double sOut        = s(Tsol, muBsol, muQsol, muSsol);
    double BOut        = B(Tsol, muBsol, muQsol, muSsol);
    double SOut        = S(Tsol, muBsol, muQsol, muSsol);
    double QOut        = Q(Tsol, muBsol, muQsol, muSsol);
    double eOut        = (sOut*Tsol + muBsol*BOut + muQsol*QOut + muSsol*SOut)/hc - POut;
    results[0]         = eOut;
    results[1]         = BOut;
    results[2]         = SOut;
    results[3]         = QOut;
  }

  void get_sBSQ( double point[], double results[] )
  {
    // point: (T, muB, muQ, muS)
    results[0] = s( hc*point[0], hc*point[1], hc*point[2], hc*point[3] );
    results[1] = B( hc*point[0], hc*point[1], hc*point[2], hc*point[3] );
    results[2] = S( hc*point[0], hc*point[1], hc*point[2], hc*point[3] );
    results[3] = Q( hc*point[0], hc*point[1], hc*point[2], hc*point[3] );
  }

  void get_full_thermo( double point[], double results[] )
  {
    // point: (T, muB, muQ, muS)
    const double Tsol   = hc*point[0], muBsol = hc*point[1],
                 muQsol = hc*point[2], muSsol = hc*point[3];
    double POut = p(Tsol, muBsol, muQsol, muSsol);
    double sOut = s(Tsol, muBsol, muQsol, muSsol);
    double BOut = B(Tsol, muBsol, muQsol, muSsol);
    double SOut = S(Tsol, muBsol, muQsol, muSsol);
    double QOut = Q(Tsol, muBsol, muQsol, muSsol);
    double eOut = (sOut*Tsol + muBsol*BOut + muQsol*QOut + muSsol*SOut)/hc - POut;

    results[0]  = POut;
    results[1]  = sOut;
    results[2]  = BOut;
    results[3]  = SOut;
    results[4]  = QOut;
    results[5]  = eOut;
    results[6]  = 1.0/3.0;  // conformal

    results[7]  = P2B2(Tsol, muBsol, muQsol, muSsol);
    results[8]  = P2Q2(Tsol, muBsol, muQsol, muSsol);
    results[9]  = P2S2(Tsol, muBsol, muQsol, muSsol);
    results[10] = P2BQ(Tsol, muBsol, muQsol, muSsol);
    results[11] = P2BS(Tsol, muBsol, muQsol, muSsol);
    results[12] = P2QS(Tsol, muBsol, muQsol, muSsol);
    results[13] = P2TB(Tsol, muBsol, muQsol, muSsol);
    results[14] = P2TQ(Tsol, muBsol, muQsol, muSsol);
    results[15] = P2TS(Tsol, muBsol, muQsol, muSsol);
    results[16] = P2T2(Tsol, muBsol, muQsol, muSsol);
  }
  
}

#endif