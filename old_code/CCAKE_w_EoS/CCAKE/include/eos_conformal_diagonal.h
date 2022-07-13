#ifndef EOS_CONFORMAL_DIAGONAL_H
#define EOS_CONFORMAL_DIAGONAL_H

#include <cmath>
#include <vector>

#include "constants.h"
#include "eos_base.h"
#include "eos_header.h"

using std::abs;
using std::pow;

class EoS_conformal_diagonal: public EoS_base
{
private:

  double c, T0, muB0, muS0, muQ0;
  static constexpr double four_thirds = 4.0/3.0;
  static constexpr double two_to_two_thirds = pow(2.0, 2.0/3.0);

  inline double sgn(double val) { return (0.0 < val) - (val < 0.0); }

public:
  // default constructor/destructor
  EoS_conformal_diagonal(){std::cout << "In default constructor!\n"; abort();}
  virtual ~EoS_conformal_diagonal(){}

  EoS_conformal_diagonal( const double c_in,
                          const double T0_in, const double muB0_in,
                          const double muS0_in, const double muQ0_in,
                          const std::vector<double> & tbqs_minima_in,
                          const std::vector<double> & tbqs_maxima_in,
                          const std::string & name_in = "conformal_diagonal")
    { c = c_in; T0 = T0_in; muB0 = muB0_in; muS0 = muS0_in; muQ0 = muQ0_in;
      tbqs_minima = tbqs_minima_in; tbqs_maxima = tbqs_maxima_in; name = name_in; }

  // this is the criterion for a guaranteed solution with this EoS
  bool eBSQ_has_solution( const double e0, const double rhoB0,
                          const double rhoS0, const double rhoQ0 )
  {
    return e0 >= 3.0*( pow(muB0*abs(rhoB0), four_thirds)
                      + pow(muS0*abs(rhoS0), four_thirds)
                      + pow(muQ0*abs(rhoQ0), four_thirds) )
                  / ( 4.0 * two_to_two_thirds * pow(c, 1.0/3.0) );
  }

  // makes an educated guess where the root should be
  std::vector<double> get_tbqs_seed_from_eBSQ(
                        const double e0, const double rhoB0,
                        const double rhoS0, const double rhoQ0 )
  {
    double den = two_to_two_thirds * pow(c, 1.0/3.0);
    double disc = e0 - 3.0*( pow(muB0*abs(rhoB0), four_thirds)
                      + pow(muS0*abs(rhoS0), four_thirds)
                      + pow(muQ0*abs(rhoQ0), four_thirds) )
                  / ( 4.0 * two_to_two_thirds * pow(c, 1.0/3.0) );
    return std::vector<double>({ T0*pow(disc/(3.0*c), 0.25),
                                sgn(rhoB0)*pow( muB0*muB0*muB0*muB0*abs(rhoB0), 1.0/3.0 )/den,
                                sgn(rhoQ0)*pow( muQ0*muQ0*muQ0*muQ0*abs(rhoQ0), 1.0/3.0 )/den,
                                sgn(rhoS0)*pow( muS0*muS0*muS0*muS0*abs(rhoS0), 1.0/3.0 )/den
                               });
  }

  // makes an educated guess where the root should be
  std::vector<double> get_tbqs_seed_from_sBSQ(
                        const double s0, const double rhoB0,
                        const double rhoS0, const double rhoQ0 )
  {
    double den = two_to_two_thirds * pow(c, 1.0/3.0);
    return std::vector<double>({ pow( T0*T0*T0*T0*s0, 1.0/3.0 )/den,
                                 sgn(rhoB0)*pow( muB0*muB0*muB0*muB0*abs(rhoB0), 1.0/3.0 )/den,
                                 sgn(rhoQ0)*pow( muQ0*muQ0*muQ0*muQ0*abs(rhoQ0), 1.0/3.0 )/den,
                                 sgn(rhoS0)*pow( muS0*muS0*muS0*muS0*abs(rhoS0), 1.0/3.0 )/den
                               });
  }

  double p(double T, double muB, double muQ, double muS)
	{
		double x = (T/T0)*(T/T0)*(T/T0)*(T/T0)
                + (muB/muB0)*(muB/muB0)*(muB/muB0)*(muB/muB0)
                + (muS/muS0)*(muS/muS0)*(muS/muS0)*(muS/muS0)
                + (muQ/muQ0)*(muQ/muQ0)*(muQ/muQ0)*(muQ/muQ0);
    double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
		return cp*x;
	}
	
	double s(double T, double muB, double muQ, double muS)
	{
    double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
		return 4.0*cp*T*T*T/(T0*T0*T0*T0);
	}

	double B(double T, double muB, double muQ, double muS)
	{
    double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
		return 4.0*cp*muB*muB*muB/(muB0*muB0*muB0*muB0);
	}

	double S(double T, double muB, double muQ, double muS)
	{
    double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
		return 4.0*cp*muS*muS*muS/(muS0*muS0*muS0*muS0);
	}

	double Q(double T, double muB, double muQ, double muS)
	{
    double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
		return 4.0*cp*muQ*muQ*muQ/(muQ0*muQ0*muQ0*muQ0);
	}

	double P2B2(double T, double muB, double muQ, double muS)
	{
    double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
		return 12.0*cp*muB*muB/(muB0*muB0*muB0*muB0);
	}
	double P2Q2(double T, double muB, double muQ, double muS)
	{
    double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
		return 12.0*cp*muQ*muQ/(muQ0*muQ0*muQ0*muQ0);
	}
	double P2S2(double T, double muB, double muQ, double muS)
	{
    double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
		return 12.0*cp*muS*muS/(muS0*muS0*muS0*muS0);
	}
	
	double P2BQ(double T, double muB, double muQ, double muS) { return TINY; }
	double P2BS(double T, double muB, double muQ, double muS) { return TINY; }
	double P2QS(double T, double muB, double muQ, double muS) { return TINY; }
	
	double P2TB(double T, double muB, double muQ, double muS) { return TINY; }
	double P2TQ(double T, double muB, double muQ, double muS) { return TINY; }
	double P2TS(double T, double muB, double muQ, double muS) { return TINY; }

	double P2T2(double T, double muB, double muQ, double muS)
	{
    double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
		return 12.0*cp*T*T/(T0*T0*T0*T0);
	}

  void get_eBSQ( double point[], double results[] )
  {
    // point: (T, muB, muQ, muS)
    const double Tsol  = point[0],
      muBsol           = point[1],
      muQsol           = point[2],
      muSsol           = point[3];
//std::cout << "LINE: " << __LINE__ << "   " << point[0] << "\n";
//std::cout << "LINE: " << __LINE__ << "   " << point[1] << "\n";
//std::cout << "LINE: " << __LINE__ << "   " << point[2] << "\n";
//std::cout << "LINE: " << __LINE__ << "   " << point[3] << "\n";
    double POut        = p(Tsol, muBsol, muQsol, muSsol);
    double sOut        = s(Tsol, muBsol, muQsol, muSsol);
    double BOut        = B(Tsol, muBsol, muQsol, muSsol);
    double SOut        = S(Tsol, muBsol, muQsol, muSsol);
    double QOut        = Q(Tsol, muBsol, muQsol, muSsol);
    double eOut        = (sOut*Tsol + muBsol*BOut + muQsol*QOut + muSsol*SOut) - POut;
    results[0]         = eOut;
    results[1]         = BOut;
    results[2]         = SOut;
    results[3]         = QOut;
//std::cout << "LINE: " << __LINE__ << "\n";
  }

  void get_sBSQ( double point[], double results[] )
  {
    // point: (T, muB, muQ, muS)
    results[0] = s( point[0], point[1], point[2], point[3] );
    results[1] = B( point[0], point[1], point[2], point[3] );
    results[2] = S( point[0], point[1], point[2], point[3] );
    results[3] = Q( point[0], point[1], point[2], point[3] );
  }

  void get_full_thermo( double point[], double results[] )
  {
    // point: (T, muB, muQ, muS)
    const double Tsol   = point[0], muBsol = point[1],
                 muQsol = point[2], muSsol = point[3];
    double POut = p(Tsol, muBsol, muQsol, muSsol);
    double sOut = s(Tsol, muBsol, muQsol, muSsol);
    double BOut = B(Tsol, muBsol, muQsol, muSsol);
    double SOut = S(Tsol, muBsol, muQsol, muSsol);
    double QOut = Q(Tsol, muBsol, muQsol, muSsol);
    double eOut = (sOut*Tsol + muBsol*BOut + muQsol*QOut + muSsol*SOut) - POut;

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