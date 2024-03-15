#ifndef EOS_CONFORMAL_COSH_H
#define EOS_CONFORMAL_COSH_H

#include <vector>

#include "constants.h"
#include "eos_base.h"
#include "eos_header.h"

class EoS_conformal_cosh: public EoS_base
{
private:
	//const double hc = constants::hbarc_MeVfm;

  double c, T0, muB0, muS0, muQ0;

public:
  // default constructor/destructor
  EoS_conformal_cosh(){}
  virtual ~EoS_conformal_cosh(){}

  EoS_conformal_cosh( const double c_in, const double T0_in, const double muB0_in,
                 const double muS0_in, const double muQ0_in,
                 const std::vector<double> & tbqs_minima_in,
                 const std::vector<double> & tbqs_maxima_in,
                 const std::string & name_in = "cosh_conformal")
    { c = c_in; T0 = T0_in; muB0 = muB0_in; muS0 = muS0_in; muQ0 = muQ0_in;
      tbqs_minima = tbqs_minima_in; tbqs_maxima = tbqs_maxima_in; 
      tbqs_minima_no_ext = tbqs_minima_in; tbqs_maxima_no_ext = tbqs_maxima_in;
      name = name_in; }


  double p(double T, double muB, double muQ, double muS)
	{
    double muEff = (muB+muQ+muS)/T;
    return c*T*T*T*T*std::cosh(muEff); 
	}
	
	double s(double T, double muB, double muQ, double muS)
	{
		double muEff = (muB+muQ+muS)/T;
    return c*T*T*T*(4*std::cosh(muEff)-muEff*std::sinh(muEff));
	}

	double B(double T, double muB, double muQ, double muS)
	{
		double muEff = (muB+muQ+muS)/T;
    return c*T*T*T*std::sinh(muEff);
	}

	double S(double T, double muB, double muQ, double muS)
	{
		double muEff = (muB+muQ+muS)/T;
    return c*T*T*T*std::sinh(muEff);
	}

	double Q(double T, double muB, double muQ, double muS)
	{
		double muEff = (muB+muQ+muS)/T;
    return c*T*T*T*std::sinh(muEff);
	}

	double P2B2(double T, double muB, double muQ, double muS)
	{
    double muEff = (muB+muQ+muS)/T;
    return c*T*T*std::cosh(muEff);
	}

	double P2Q2(double T, double muB, double muQ, double muS)
	{
    double muEff = (muB+muQ+muS)/T;
    return c*T*T*std::cosh(muEff);
	}
	double P2S2(double T, double muB, double muQ, double muS)
	{
    double muEff = (muB+muQ+muS)/T;
    return c*T*T*std::cosh(muEff);
	}
	
	double P2BQ(double T, double muB, double muQ, double muS)
      {  double muEff = (muB+muQ+muS)/T;
    return c*T*T*std::cosh(muEff);}
	double P2BS(double T, double muB, double muQ, double muS)
			{     double muEff = (muB+muQ+muS)/T;
    return c*T*T*std::cosh(muEff); }
	double P2QS(double T, double muB, double muQ, double muS)
			{     double muEff = (muB+muQ+muS)/T;
    return c*T*T*std::cosh(muEff); }
	
	double P2TB(double T, double muB, double muQ, double muS)
			{ 
        double muEff = (muB+muQ+muS)/T;
        return c*T*T*(3.*std::sinh(muEff)-muEff*std::cosh(muEff));
      }
	double P2TQ(double T, double muB, double muQ, double muS)
			{ double muEff = (muB+muQ+muS)/T;
        return c*T*T*(3.*std::sinh(muEff)-muEff*std::cosh(muEff)); }
	double P2TS(double T, double muB, double muQ, double muS)
			{ double muEff = (muB+muQ+muS)/T;
        return c*T*T*(3.*std::sinh(muEff)-muEff*std::cosh(muEff)); }
	double P2T2(double T, double muB, double muQ, double muS)
	{
    double muEff = (muB+muQ+muS)/T;
    return c*T*T*(muEff*muEff*std::cosh(muEff)-6.*muEff*std::sinh(muEff)+12.*std::cosh(muEff));
	}

  void get_eBSQ( double point[], double results[] )
  {
    // point: (T, muB, muQ, muS)
    const double Tsol   = point[0], muBsol = point[1],
                 muQsol = point[2], muSsol = point[3];
//std::cout << __FUNCTION__ << "::" << __LINE__ << ": " << point[0] << std::endl;
//std::cout << __FUNCTION__ << "::" << __LINE__ << ": " << point[1] << std::endl;
//std::cout << __FUNCTION__ << "::" << __LINE__ << ": " << point[2] << std::endl;
//std::cout << __FUNCTION__ << "::" << __LINE__ << ": " << point[3] << std::endl;
//std::cout << __FUNCTION__ << "::" << __LINE__ << ": " << results[0] << std::endl;
//std::cout << __FUNCTION__ << "::" << __LINE__ << ": " << results[1] << std::endl;
//std::cout << __FUNCTION__ << "::" << __LINE__ << ": " << results[2] << std::endl;
//std::cout << __FUNCTION__ << "::" << __LINE__ << ": " << results[3] << std::endl;
    double POut        = p(Tsol, muBsol, muQsol, muSsol);
//std::cout << __FUNCTION__ << "::" << __LINE__ << std::endl;
    double sOut        = s(Tsol, muBsol, muQsol, muSsol);
//std::cout << __FUNCTION__ << "::" << __LINE__ << std::endl;
    double BOut        = B(Tsol, muBsol, muQsol, muSsol);
//std::cout << __FUNCTION__ << "::" << __LINE__ << std::endl;
    double SOut        = S(Tsol, muBsol, muQsol, muSsol);
//std::cout << __FUNCTION__ << "::" << __LINE__ << std::endl;
    double QOut        = Q(Tsol, muBsol, muQsol, muSsol);
//std::cout << __FUNCTION__ << "::" << __LINE__ << std::endl;
    double eOut        = (sOut*Tsol + muBsol*BOut + muQsol*QOut + muSsol*SOut) - POut;
//std::cout << __FUNCTION__ << "::" << __LINE__ << std::endl;
    results[0]         = eOut;
    results[1]         = BOut;
    results[2]         = SOut;
    results[3]         = QOut;
//std::cout << __FUNCTION__ << "::" << __LINE__ << std::endl;
    return;
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