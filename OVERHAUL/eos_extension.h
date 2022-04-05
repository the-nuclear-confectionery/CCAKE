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

  double b0000, b0002, b0020, b0200, b2000;
  double b0022, b0202, b2002, b0220, b2020, b2200;
  double b0004, b0040, b0400, b4000;

  //////////////////////////////////////////////////////////////////////////////
  // function definitions
  void set_coeffs( const double point[], double thermodynamics[] );

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

    b2000 = (3.0*s0 - T0*chiTT0 - muB0*chiTB0 - muQ0*chiTQ0 - muS0*chiTS0) / (4.0*T0+TINY);
    b0200 = (3.0*rhoB0 - T0*chiTB0 - muB0*chiBB0 - muQ0*chiBQ0 - muS0*chiBS0) / (4.0*muB0+TINY);
    b0020 = (3.0*rhoS0 - T0*chiTS0 - muB0*chiBS0 - muQ0*chiSQ0 - muS0*chiSS0) / (4.0*muS0+TINY);
    b0002 = (3.0*rhoQ0 - T0*chiTQ0 - muB0*chiBQ0 - muQ0*chiQQ0 - muS0*chiSQ0) / (4.0*muQ0+TINY);

    b2200 = chiTB0 / (4.0*T0*muB0+TINY);
    b2020 = chiTS0 / (4.0*T0*muS0+TINY);
    b2002 = chiTQ0 / (4.0*T0*muQ0+TINY);
    b0220 = chiBS0 / (4.0*muB0*muS0+TINY);
    b0202 = chiBQ0 / (4.0*muB0*muQ0+TINY);
    b0022 = chiSQ0 / (4.0*muS0*muQ0+TINY);

    b4000 = (T0*chiTT0 - s0) / (8.0*T0*T0*T0+TINY);
    b0400 = (muB0*chiBB0 - rhoB0) / (8.0*muB0*muB0*muB0+TINY);
    b0040 = (muS0*chiSS0 - rhoS0) / (8.0*muS0*muS0*muS0+TINY);
    b0004 = (muQ0*chiQQ0 - rhoQ0) / (8.0*muQ0*muQ0*muQ0+TINY);

cout << __PRETTY_FUNCTION__ << "::" << __LINE__ << endl;
cout << T0 << endl;
cout << muB0 << endl;
cout << muQ0 << endl;
cout << muS0 << endl;
cout << "------------------" << endl;
cout << b0000 << endl;
cout << "------------------" << endl;
cout << b2000 << endl;
cout << b0200 << endl;
cout << b0020 << endl;
cout << b0002 << endl;
cout << "------------------" << endl;
cout << b4000 << endl;
cout << b0400 << endl;
cout << b0040 << endl;
cout << b0004 << endl;
cout << "------------------" << endl;
cout << b2200 << endl;
cout << b2020 << endl;
cout << b2002 << endl;
cout << b0220 << endl;
cout << b0202 << endl;
cout << b0022 << endl;
  }

  double p(double T, double muB, double muQ, double muS);
cout << __PRETTY_FUNCTION__ << "::" << __LINE__ << endl;
cout << T << endl;
cout << muB << endl;
cout << muQ << endl;
cout << muS << endl;
cout << "------------------" << endl;
cout << b0000 << endl;
cout << "------------------" << endl;
cout << b2000 << endl;
cout << b0200 << endl;
cout << b0020 << endl;
cout << b0002 << endl;
cout << "------------------" << endl;
cout << b4000 << endl;
cout << b0400 << endl;
cout << b0040 << endl;
cout << b0004 << endl;
cout << "------------------" << endl;
cout << b2200 << endl;
cout << b2020 << endl;
cout << b2002 << endl;
cout << b0220 << endl;
cout << b0202 << endl;
cout << b0022 << endl;
		return b0000 + b2000*T*T + b0200*muB*muB + b0020*muS*muS + b0002*muQ*muQ
            + b4000*T*T*T*T + b0400*muB*muB*muB*muB
            + b0040*muS*muS*muS*muS + b0004*muQ*muQ*muQ*muQ
            + b2200*T*T*muB*muB + b2020*T*T*muS*muS + b2002*T*T*muQ*muQ
            + b0220*muB*muB*muS*muS + b0202*muB*muB*muQ*muQ + b0022*muS*muS*muQ*muQ;
	}
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