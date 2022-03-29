#ifndef EOS_EXTENSION_H
#define EOS_EXTENSION_H

#include <algorithm>
#include <cmath>

#include "constants.h"

namespace eos_extension
{
  using namespace constants;

	const double hc = hbarc_MeVfm;

  double b0000, b0002, b0020, b0200, b2000;
  double b0022, b0202, b2002, b0220, b2020, b2200;
  double b0004, b0040, b0400, b4000;

  //////////////////////////////////////////////////////////////////////////////
  // function prototypes
  void get_nonconformal_extension( double point[], double thermodynamics[] );
  void set_coeffs( double point[], double thermodynamics[] );
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
  void get_full_thermo( double point[], double results[] );
  void project_to_boundary( double point[], double minima[], double maxima[] );

  //////////////////////////////////////////////////////////////////////////////
  // function definitions
  void get_nonconformal_extension( double point[], double thermodynamics[] )
  {
    // determine parameters from thermodynamic quantities
    set_coeffs( point, thermodynamics );

    // evaluate extension and return result (stored in thermodynamics)
    get_full_thermo( point, thermodynamics );
  }

  void set_coeffs( double point[], double thermodynamics[] )
  {

    double T0 = point[0], muB0 = point[1], muQ0 = point[2], muS0 = point[3];

    double p0     = thermodynamics[0];
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

    b0000 = p0 - (5.0/8.0)*(T0*s0+muB0*rhoB0+muS0*rhoS0+muQ0*rhoQ0)
            + ( chiTT0*T0*T0 + chiTB0*T0*rhoB0 + chiTS0*T0*rhoS0 + chiTQ0*T0*rhoQ0
                + chiBB0*rhoB0*rhoB0 + chiBS0*rhoB0*rhoS0 + chiBQ0*rhoB0*rhoQ0
                + chiSS0*rhoS0*rhoS0 + chiSQ0*rhoS0*rhoQ0 + chiQQ0*rhoQ0*rhoQ0 ) / 8.0;

    b2000 = (3.0*s0 - T0*chiTT0 - muB0*chiTB0 - muQ0*chiTQ0 - muS0*chiTS0) / (4.0*T0);
    b0200 = (3.0*rhoB0 - T0*chiTB0 - muB0*chiBB0 - muQ0*chiBQ0 - muS0*chiBS0) / (4.0*muB0);
    b0020 = (3.0*rhoS0 - T0*chiTS0 - muB0*chiBS0 - muQ0*chiSQ0 - muS0*chiSS0) / (4.0*muS0);
    b0002 = (3.0*rhoQ0 - T0*chiTQ0 - muB0*chiBQ0 - muQ0*chiQQ0 - muS0*chiSQ0) / (4.0*muQ0);

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
		return b0000 + b2000*T*T + b0200*muB*muB + b0020*muS*muS + b0002*muQ*muQ
            + b4000*T*T*T*T + b0400*muB*muB*muB*muB
            + b0040*muS*muS*muS*muS + b0004*muQ*muQ*muQ*muQ
            + b2200*T*T*muB*muB + b2020*T*T*muS*muS + b2002*T*T*muQ*muQ
            + b0220*muB*muB*muS*muS + b0202*muB*muB*muQ*muQ + b0022*muS*muS*muQ*muQ;
	}
	
	double s(double T, double muB, double muQ, double muS)
	{
		2.0*T*(b2000 + 2.0*T*T*b4000 + muB*muB*b2200 + muQ*muQ*b2002 + muS*muS*b2020);
	}

	double B(double T, double muB, double muQ, double muS)
	{
		2.0*muB*(b0200 + T*T*b2200 + 2.0*muB*muB*b0400 + muS*muS*b0220 + muQ*muQ*b0202);
	}

	double S(double T, double muB, double muQ, double muS)
	{
		2.0*T*(b0020 + T*T*b2020 + muB*muB*b0220 + 2.0*muS*muS*b0040 + muQ*muQ*b0022);
	}

	double Q(double T, double muB, double muQ, double muS)
	{
		2.0*T*(b0002 + T*T*b2002 + muB*muB*b0202 + muS*muS*b0022 + 2.0*muQ*muQ*b0004);
	}

	double P2T2(double T, double muB, double muQ, double muS)
	{
		2.0*(b2000 + 6.0*T*T*b4000 + muB*muB*b2200 + muS*muS*b2020 + muQ*muQ*b2002);
	}
	double P2B2(double T, double muB, double muQ, double muS)
	{
		2.0*(b0200 + T*T*b2200 + 6.0*muB*muB*b0400 + muS*muS*b0220 + muQ*muQ*b0202);
	}
	double P2S2(double T, double muB, double muQ, double muS)
	{
		2.0*(b0020 + T*T*b2020 + muB*muB*b0220 + 6.0*muS*muS*b0040 + muQ*muQ*b0022);
	}
	double P2Q2(double T, double muB, double muQ, double muS)
	{
		2.0*(b0002 + T*T*b2002 + muB*muB*b0202 + muS*muS*b0022 + 6.0*muQ*muQ*b0004);
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
    //results[6]  = 1.0/3.0;  // conformal

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

    ////////////////////////////////////////////////////////////////////////////

    const double T = point[0], muB = point[1], muQ = point[2], muS = point[3];

    double C1T = results[1];
    double C1B = results[2];
    double C1Q = results[3];
    double C1S = results[4];

    double C2B2 = results[7];
    double C2Q2 = results[8];
    double C2S2 = results[9];
    double C2BQ = results[10];
    double C2BS = results[11];
    double C2QS = results[12];
    double C2TB = results[13];
    double C2TQ = results[14];
    double C2TS = results[15];
    double C2T2 = results[16];

    // speed of sound
    results[6] = T*(-(C2BQ*C2S2*C2TQ*C1B) - C2BQ*C2S2*C2TB*C1Q - pow(C2BS,2)*C2TQ*C1Q + C2B2*C2S2*C2TQ*C1Q + C2BQ*C2BS*C2TS*C1Q + C2BQ*C2BS*C2TQ*C1S - pow(C2BQ,2)*C2TS*C1S + pow(C2BQ,2)*C2S2*C1T 
              + pow(C2QS,2)*(-(C2TB*C1B) + C2B2*C1T) + C2Q2*(C2S2*C2TB*C1B - C2BS*C2TS*C1B - C2BS*C2TB*C1S + C2B2*C2TS*C1S + pow(C2BS,2)*C1T - C2B2*C2S2*C1T) + C2QS*(C2BQ*C2TS*C1B - C2B2*C2TS*C1Q 
              + C2BQ*C2TB*C1S - C2B2*C2TQ*C1S + C2BS*(C2TQ*C1B + C2TB*C1Q - 2.0*C2BQ*C1T)))
                 *1.0/((pow(C2BQ,2)*C2S2*C2T2 - pow(C2QS,2)*pow(C2TB,2) + C2Q2*C2S2*pow(C2TB,2) - 2.0*C2BQ*C2S2*C2TB*C2TQ + pow(C2BS,2)*(C2Q2*C2T2 - pow(C2TQ,2)) 
              + 2.0*C2BQ*C2QS*C2TB*C2TS - pow(C2BQ,2)*pow(C2TS,2) + 2.0*C2BS*(-(C2BQ*C2QS*C2T2) + C2QS*C2TB*C2TQ - C2Q2*C2TB*C2TS + C2BQ*C2TQ*C2TS) + C2B2*(pow(C2QS,2)*C2T2 - C2Q2*C2S2*C2T2 
              + C2S2*pow(C2TQ,2) - 2.0*C2QS*C2TQ*C2TS + C2Q2*pow(C2TS,2)))*T) 
            + T*(C1B/(muB*C1B + muQ*C1Q + muS*C1S + T*C1T))*(-(C2QS*muB*(C2BQ*(C2TS*C1B + C2TB*C1S) - C2B2*(C2TS*C1Q + C2TQ*C1S) + C2BS*(C2TQ*C1B + C2TB*C1Q - 2*C2BQ*C1T))) + muB*((pow(C2BS,2) 
              - C2B2*C2S2)*C2TQ*C1Q + C2BQ*(C2S2*C2TQ*C1B + C2S2*C2TB*C1Q - C2BS*C2TS*C1Q - C2BS*C2TQ*C1S) + pow(C2BQ,2)*(C2TS*C1S - C2S2*C1T) + C2Q2*(-(C2S2*C2TB*C1B) + C2BS*C2TS*C1B 
              + C2BS*C2TB*C1S - C2B2*C2TS*C1S - pow(C2BS,2)*C1T + C2B2*C2S2*C1T)) + C2QS*(-2*C2TQ*C2TS*C1B + C2TB*C2TS*C1Q + C2TB*C2TQ*C1S - C2T2*(C2BS*C1Q + C2BQ*C1S) + C2BS*C2TQ*C1T 
              + C2BQ*C2TS*C1T)*T + (-((C2BS*C2TQ - C2BQ*C2TS)*(-(C2TS*C1Q) + C2TQ*C1S)) + C2S2*(C2BQ*C2T2*C1Q + C2TQ*(C2TQ*C1B - C2TB*C1Q - C2BQ*C1T)) + C2Q2*(C2BS*C2T2*C1S + C2TS*(C2TS*C1B 
              - C2TB*C1S - C2BS*C1T) + C2S2*(-(C2T2*C1B) + C2TB*C1T)))*T + pow(C2QS,2)*(C2TB*muB*C1B - C2B2*muB*C1T + C2T2*C1B*T - C2TB*C1T*T))/((pow(C2BQ,2)*C2S2*C2T2 - pow(C2QS,2)*pow(C2TB,2) 
              + C2Q2*C2S2*pow(C2TB,2) - 2*C2BQ*C2S2*C2TB*C2TQ + pow(C2BS,2)*(C2Q2*C2T2 - pow(C2TQ,2)) + 2*C2BQ*C2QS*C2TB*C2TS - pow(C2BQ,2)*pow(C2TS,2) + 2*C2BS*(-(C2BQ*C2QS*C2T2) + C2QS*C2TB*C2TQ 
              - C2Q2*C2TB*C2TS + C2BQ*C2TQ*C2TS) + C2B2*(pow(C2QS,2)*C2T2 - C2Q2*C2S2*C2T2 + C2S2*pow(C2TQ,2) - 2*C2QS*C2TQ*C2TS + C2Q2*pow(C2TS,2)))*T)
            + T*(C1Q/(muB*C1B + muQ*C1Q + muS*C1S + T*C1T))*(pow(C2QS,2)*muQ*(C2TB*C1B - C2B2*C1T) - C2QS*muQ*(C2BQ*(C2TS*C1B + C2TB*C1S) - C2B2*(C2TS*C1Q + C2TQ*C1S) + C2BS*(C2TQ*C1B + C2TB*C1Q 
              - 2*C2BQ*C1T)) + muQ*((pow(C2BS,2) - C2B2*C2S2)*C2TQ*C1Q + C2BQ*(C2S2*C2TQ*C1B + C2S2*C2TB*C1Q - C2BS*C2TS*C1Q - C2BS*C2TQ*C1S) + pow(C2BQ,2)*(C2TS*C1S - C2S2*C1T) + C2Q2*(-(C2S2*C2TB*C1B) 
              + C2BS*C2TS*C1B + C2BS*C2TB*C1S - C2B2*C2TS*C1S - pow(C2BS,2)*C1T + C2B2*C2S2*C1T)) + C2QS*(-(C2BS*C2T2*C1B) + C2TB*C2TS*C1B + C2B2*C2T2*C1S - pow(C2TB,2)*C1S + C2BS*C2TB*C1T 
              - C2B2*C2TS*C1T)*T + (C2BS*C2TQ*C2TS*C1B + pow(C2BS,2)*C2T2*C1Q - 2*C2BS*C2TB*C2TS*C1Q + C2B2*pow(C2TS,2)*C1Q + C2BS*C2TB*C2TQ*C1S - C2B2*C2TQ*C2TS*C1S - pow(C2BS,2)*C2TQ*C1T 
              + C2S2*(-(C2TB*C2TQ*C1B) - C2B2*C2T2*C1Q + pow(C2TB,2)*C1Q + C2B2*C2TQ*C1T) + C2BQ*(-(C2BS*C2T2*C1S) + C2TS*(-(C2TS*C1B) + C2TB*C1S + C2BS*C1T) + C2S2*(C2T2*C1B - C2TB*C1T)))*T)
                 *1.0/((pow(C2BQ,2)*C2S2*C2T2 - pow(C2QS,2)*pow(C2TB,2) + C2Q2*C2S2*pow(C2TB,2) - 2*C2BQ*C2S2*C2TB*C2TQ + pow(C2BS,2)*(C2Q2*C2T2 - pow(C2TQ,2)) + 2*C2BQ*C2QS*C2TB*C2TS 
              - pow(C2BQ,2)*pow(C2TS,2) + 2*C2BS*(-(C2BQ*C2QS*C2T2) + C2QS*C2TB*C2TQ - C2Q2*C2TB*C2TS + C2BQ*C2TQ*C2TS) + C2B2*(pow(C2QS,2)*C2T2 - C2Q2*C2S2*C2T2 + C2S2*pow(C2TQ,2) 
              - 2*C2QS*C2TQ*C2TS + C2Q2*pow(C2TS,2)))*T)
            + T*(C1S/(muB*C1B + muQ*C1Q + muS*C1S + T*C1T))*(pow(C2QS,2)*muS*(C2TB*C1B - C2B2*C1T) - C2QS*muS*(C2BQ*(C2TS*C1B + C2TB*C1S) - C2B2*(C2TS*C1Q + C2TQ*C1S) + C2BS*(C2TQ*C1B + C2TB*C1Q 
              - 2*C2BQ*C1T)) + muS*((pow(C2BS,2) - C2B2*C2S2)*C2TQ*C1Q + C2BQ*(C2S2*C2TQ*C1B + C2S2*C2TB*C1Q - C2BS*C2TS*C1Q - C2BS*C2TQ*C1S) + pow(C2BQ,2)*(C2TS*C1S - C2S2*C1T) 
              + C2Q2*(-(C2S2*C2TB*C1B) + C2BS*C2TS*C1B + C2BS*C2TB*C1S - C2B2*C2TS*C1S - pow(C2BS,2)*C1T + C2B2*C2S2*C1T)) + C2QS*(-(C2BQ*C2T2*C1B) + C2TB*C2TQ*C1B + C2B2*C2T2*C1Q - pow(C2TB,2)*C1Q 
              + C2BQ*C2TB*C1T - C2B2*C2TQ*C1T)*T + (C2BQ*C2TQ*C2TS*C1B + C2BQ*C2TB*C2TS*C1Q - C2B2*C2TQ*C2TS*C1Q + pow(C2BQ,2)*C2T2*C1S - 2*C2BQ*C2TB*C2TQ*C1S + C2B2*pow(C2TQ,2)*C1S 
              - pow(C2BQ,2)*C2TS*C1T + C2Q2*(-(C2TB*C2TS*C1B) - C2B2*C2T2*C1S + pow(C2TB,2)*C1S + C2B2*C2TS*C1T) + C2BS*(-(C2BQ*C2T2*C1Q) + C2TQ*(-(C2TQ*C1B) + C2TB*C1Q + C2BQ*C1T) 
              + C2Q2*(C2T2*C1B - C2TB*C1T)))*T)/((pow(C2BQ,2)*C2S2*C2T2 - pow(C2QS,2)*pow(C2TB,2) + C2Q2*C2S2*pow(C2TB,2) - 2*C2BQ*C2S2*C2TB*C2TQ + pow(C2BS,2)*(C2Q2*C2T2 - pow(C2TQ,2)) 
              + 2*C2BQ*C2QS*C2TB*C2TS - pow(C2BQ,2)*pow(C2TS,2) + 2*C2BS*(-(C2BQ*C2QS*C2T2) + C2QS*C2TB*C2TQ - C2Q2*C2TB*C2TS + C2BQ*C2TQ*C2TS) + C2B2*(pow(C2QS,2)*C2T2 - C2Q2*C2S2*C2T2 
              + C2S2*pow(C2TQ,2) - 2*C2QS*C2TQ*C2TS + C2Q2*pow(C2TS,2)))*T);

    ////////////////////////////////////////////////////////////////////////////


  }

  void project_to_boundary( double point[], double minima[], double maxima[] )
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