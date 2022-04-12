#ifndef EOS_CONFORMAL_DIAGONAL_H
#define EOS_CONFORMAL_DIAGONAL_H

#include <vector>

#include "constants.h"
#include "eos_base.h"
#include "eos_header.h"

class EoS_conformal_diagonal: public EoS_base
{
private:
	const double hc = constants::hbarc_MeVfm;

  double c, T0, muB0, muS0, muQ0;

public:
  // default constructor/destructor
  EoS_conformal_diagonal(){}
  virtual ~EoS_conformal_diagonal(){}

  EoS_conformal_diagonal( const double c_in,
                          const double T0_in, const double muB0_in,
                          const double muS0_in, const double muQ0_in,
                          const std::vector<double> & tbqs_minima_in,
                          const std::vector<double> & tbqs_maxima_in,
                          const std::string & name_in = "conformal_diagonal")
    { c = c_in; T0 = T0_in; muB0 = muB0_in; muS0 = muS0_in; muQ0 = muQ0_in;
      tbqs_minima = tbqs_minima_in; tbqs_maxima = tbqs_maxima_in; name = name_in; }


  double p(double T, double muB, double muQ, double muS)
	{
		double x = (T/T0)*(T/T0)*(T/T0)*(T/T0)
                + (muB/muB0)*(muB/muB0)*(muB/muB0)*(muB/muB0)
                + (muS/muS0)*(muS/muS0)*(muS/muS0)*(muS/muS0)
                + (muQ/muQ0)*(muQ/muQ0)*(muQ/muQ0)*(muQ/muQ0);
    double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
		return cp*x;
	}
	
	double s(double T, double muB, double muQ, double muS)
	{
    double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
		return 4.0*cp*hc*T*T*T/(T0*T0*T0*T0);
	}

	double B(double T, double muB, double muQ, double muS)
	{
    double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
		return 4.0*cp*hc*muB*muB*muB/(muB0*muB0*muB0*muB0);
	}

	double S(double T, double muB, double muQ, double muS)
	{
    double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
		return 4.0*cp*hc*muS*muS*muS/(muS0*muS0*muS0*muS0);
	}

	double Q(double T, double muB, double muQ, double muS)
	{
    double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
		return 4.0*cp*hc*muQ*muQ*muQ/(muQ0*muQ0*muQ0*muQ0);
	}

	double P2B2(double T, double muB, double muQ, double muS)
	{
    double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
		return 12.0*cp*hc*hc*muB*muB/(muB0*muB0*muB0*muB0);
	}
	double P2Q2(double T, double muB, double muQ, double muS)
	{
    double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
		return 12.0*cp*hc*hc*muQ*muQ/(muQ0*muQ0*muQ0*muQ0);
	}
	double P2S2(double T, double muB, double muQ, double muS)
	{
    double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
		return 12.0*cp*hc*hc*muS*muS/(muS0*muS0*muS0*muS0);
	}
	
	double P2BQ(double T, double muB, double muQ, double muS) { return 0.0; }
	double P2BS(double T, double muB, double muQ, double muS) { return 0.0; }
	double P2QS(double T, double muB, double muQ, double muS) { return 0.0; }
	
	double P2TB(double T, double muB, double muQ, double muS) { return 0.0; }
	double P2TQ(double T, double muB, double muQ, double muS) { return 0.0; }
	double P2TS(double T, double muB, double muQ, double muS) { return 0.0; }

	double P2T2(double T, double muB, double muQ, double muS)
	{
    double cp = c*T0*T0*T0*T0/(hc*hc*hc*hc);  // c dimensionless, cp - 1/fm^4
		return 12.0*hc*hc*cp*T*T/(T0*T0*T0*T0);
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
  
};

#endif