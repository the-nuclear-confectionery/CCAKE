#ifndef EOS_CONFORMAL_EXTENSION_H
#define EOS_CONFORMAL_EXTENSION_H

#include <algorithm>
#include <cmath>

#include "constants.h"

namespace eos_conformal_extension
{
  using namespace constants;

	const double hc = hbarc_MeVfm;

  //double b0000, b0002, b0020, b0200, b2000;
  double b0022, b0202, b2002, b0220, b2020, b2200;
  double b0004, b0040, b0400, b4000;

  //////////////////////////////////////////////////////////////////////////////
  // function definitions
  void set_coeffs( double point[], double thermodynamics[] )
  {

    double T0 = point[0], muB0 = point[1], muQ0 = point[2], muS0 = point[3];

    //double p0     = thermodynamics[0];
    double s0     = thermodynamics[1];
    double rhoB0  = thermodynamics[2];
    double rhoS0  = thermodynamics[3];
    double rhoQ0  = thermodynamics[4];
    //double e0     = thermodynamics[5];
    //double cs20   = thermodynamics[6];
    double chiBB0 = thermodynamics[7];
    double chiQQ0 = thermodynamics[8];
    double chiSS0 = thermodynamics[9];
    double chiBQ0 = thermodynamics[10];
    double chiBS0 = thermodynamics[11];
    double chiSQ0 = thermodynamics[12];
    double chiTB0 = thermodynamics[13];
    double chiTQ0 = thermodynamics[14];
    double chiTS0 = thermodynamics[15];
    double chiTT0 = thermodynamics[16];

    b2200 = chiTB0 / (4.0*T0*muB0);
    b2020 = chiTS0 / (4.0*T0*muS0);
    b2002 = chiTQ0 / (4.0*T0*muQ0);
    b0220 = chiBS0 / (4.0*muB0*muS0);
    b0202 = chiBQ0 / (4.0*muB0*muQ0);
    b0022 = chiSQ0 / (4.0*muS0*muQ0);

    b4000 = (T0*chiTT0 - s0) / (8.0*T0*T0*T0);
    b0400 = (muB0*chiBB0 - rhoB0) / (8.0*muB0*muB0*muB0);
    b0040 = (muS0*chiSS0 - rhoS0) / (8.0*muS0*muS0*muS0);
    b0004 = (muQ0*chiQQ0 - rhoQ0) / (8.0*muQ0*muQ0*muQ0);
  }


  double p(double T, double muB, double muQ, double muS)
	{
		return b4000*T*T*T*T + b0400*muB*muB*muB*muB
            + b0040*muS*muS*muS*muS + b0004*muQ*muQ*muQ*muQ
            + b2200*T*T*muB*muB + b2020*T*T*muS*muS + b2002*T*T*muQ*muQ
            + b0220*muB*muB*muS*muS + b0202*muB*muB*muQ*muQ + b0022*muS*muS*muQ*muQ;
	}
	
	double s(double T, double muB, double muQ, double muS)
	{
		2.0*T*(2.0*T*T*b4000 + muB*muB*b2200 + muQ*muQ*b2002 + muS*muS*b2020);
	}

	double B(double T, double muB, double muQ, double muS)
	{
		2.0*muB*(T*T*b2200 + 2.0*muB*muB*b0400 + muS*muS*b0220 + muQ*muQ*b0202);
	}

	double S(double T, double muB, double muQ, double muS)
	{
		2.0*T*(T*T*b2020 + muB*muB*b0220 + 2.0*muS*muS*b0040 + muQ*muQ*b0022);
	}

	double Q(double T, double muB, double muQ, double muS)
	{
		2.0*T*(T*T*b2002 + muB*muB*b0202 + muS*muS*b0022 + 2.0*muQ*muQ*b0004);
	}

	double P2T2(double T, double muB, double muQ, double muS)
	{
		2.0*(6.0*T*T*b4000 + muB*muB*b2200 + muS*muS*b2020 + muQ*muQ*b2002);
	}
	double P2B2(double T, double muB, double muQ, double muS)
	{
		2.0*(T*T*b2200 + 6.0*muB*muB*b0400 + muS*muS*b0220 + muQ*muQ*b0202);
	}
	double P2S2(double T, double muB, double muQ, double muS)
	{
		2.0*(T*T*b2020 + muB*muB*b0220 + 6.0*muS*muS*b0040 + muQ*muQ*b0022);
	}
	double P2Q2(double T, double muB, double muQ, double muS)
	{
		2.0*(T*T*b2002 + muB*muB*b0202 + muS*muS*b0022 + 6.0*muQ*muQ*b0004);
	}
	
	double P2BQ(double T, double muB, double muQ, double muS)
      { return 4.0*muB*muQ*b0202; }
	double P2BS(double T, double muB, double muQ, double muS)
      { return 4.0*muB*muS*b0220; }
	double P2QS(double T, double muB, double muQ, double muS)
      { return 4.0*muQ*muS*b0022; }
	
	double P2TB(double T, double muB, double muQ, double muS)
      { return 4.0*T*muB*b2200; }
	double P2TS(double T, double muB, double muQ, double muS)
      { return 4.0*T*muS*b2020; }
	double P2TQ(double T, double muB, double muQ, double muS)
      { return 4.0*T*muQ*b2002; }


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
    double eOut = sOut*Tsol + muBsol*BOut + muQsol*QOut + muSsol*SOut - POut;

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

  void get_conformal_extension( double point[], double thermodynamics[] )
  {
    std::cout << "Getting conformal extension!" << std::endl;

    // determine parameters from thermodynamic quantities
    set_coeffs( point, thermodynamics );

    // evaluate extension and return result (stored in thermodynamics)
    get_full_thermo( point, thermodynamics );
  }

  void project_to_boundary( double point[], const double minima[], const double maxima[] )
  {
    double max_ratio = 0.0;
    for (int i = 0; i < 4; i++)
    {
      double extremum = point[i] < 0.0 ? minima[i] : maxima[i]; // should never be zero
      max_ratio = std::max( max_ratio, std::abs(point[i] / extremum) );
    }
    for (int i = 0; i < 4; i++) point[i] /= max_ratio;
  }

  
}

#endif