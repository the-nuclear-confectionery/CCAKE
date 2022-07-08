#include <algorithm>
#include <cmath>
#include <iostream>

#include "../include/constants.h"
#include "../include/eos_extension.h"

namespace eos_extension
{
  using namespace constants;
  using std::cout;
  using std::endl;

  double b0000 = 0.0, b0001 = 0.0, b0010 = 0.0, b0100 = 0.0, b1000 = 0.0;
  double b0011 = 0.0, b0101 = 0.0, b1001 = 0.0, b0110 = 0.0, b1010 = 0.0, b1100 = 0.0;
  double b0002 = 0.0, b0020 = 0.0, b0200 = 0.0, b2000 = 0.0;

  //////////////////////////////////////////////////////////////////////////////
  // function definitions
  void set_coeffs( const double point[], double thermodynamics[] )
  {

    double T0 = point[0], muB0 = point[1], muQ0 = point[2], muS0 = point[3];

    double p0     = thermodynamics[0];
    double s0     = thermodynamics[1];
    double rhoB0  = thermodynamics[2];
    double rhoS0  = thermodynamics[3];
    double rhoQ0  = thermodynamics[4];
    double e0     = thermodynamics[5];
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

    double p0_1_2 = sqrt(p0);
    double p0_3_2 = p0_1_2*p0_1_2*p0_1_2;

    double tr  = T0*T0*chiTT0 + muB0*muB0*chiBB0 + muS0*muS0*chiSS0 + muQ0*muQ0*chiQQ0
                 + 2.0*( T0*muB0*chiTB0 + T0*muS0*chiTS0 + T0*muQ0*chiTQ0
                        + muB0*muQ0*chiBQ0 + muB0*muS0*chiBS0 + muS0*muQ0*chiSQ0 );

    double trT = T0*chiTT0 + muB0*chiTB0 + muS0*chiTS0 + muQ0*chiTQ0;
    double trB = T0*chiTB0 + muB0*chiBB0 + muS0*chiBS0 + muQ0*chiBQ0;
    double trS = T0*chiTS0 + muB0*chiBS0 + muS0*chiSS0 + muQ0*chiSQ0;
    double trQ = T0*chiTQ0 + muB0*chiBQ0 + muS0*chiSQ0 + muQ0*chiQQ0;

    b0000 = (e0*e0 + 6.0*e0*p0 - p0*(3.0*p0 + 2.0*tr)) / (8.0*p0_3_2);

    b1000 = (2.0*p0*trT - s0   *(e0 + 3.0*p0)) / (4.0*p0_3_2);
    b0100 = (2.0*p0*trB - rhoB0*(e0 + 3.0*p0)) / (4.0*p0_3_2);
    b0010 = (2.0*p0*trS - rhoS0*(e0 + 3.0*p0)) / (4.0*p0_3_2);
    b0001 = (2.0*p0*trQ - rhoQ0*(e0 + 3.0*p0)) / (4.0*p0_3_2);

    b1100 = (s0*rhoB0    - 2.0*p0*chiTB0) / (4.0*p0_3_2);
    b1010 = (s0*rhoS0    - 2.0*p0*chiTS0) / (4.0*p0_3_2);
    b1001 = (s0*rhoQ0    - 2.0*p0*chiTQ0) / (4.0*p0_3_2);
    b0110 = (rhoB0*rhoS0 - 2.0*p0*chiBS0) / (4.0*p0_3_2);
    b0101 = (rhoB0*rhoQ0 - 2.0*p0*chiBQ0) / (4.0*p0_3_2);
    b0011 = (rhoS0*rhoQ0 - 2.0*p0*chiSQ0) / (4.0*p0_3_2);

    b2000 = (s0*s0       - 2.0*p0*chiTT0) / (8.0*p0_3_2);
    b0200 = (rhoB0*rhoB0 - 2.0*p0*chiBB0) / (8.0*p0_3_2);
    b0020 = (rhoS0*rhoS0 - 2.0*p0*chiSS0) / (8.0*p0_3_2);
    b0002 = (rhoQ0*rhoQ0 - 2.0*p0*chiQQ0) / (8.0*p0_3_2);

//std::cout << "--------------------------------------------------------------------" << std::endl;
//std::cout << "Check b0000: " << b0000 << std::endl;
//std::cout << "Check b1000: " << b1000 << std::endl;
//std::cout << "Check b0100: " << b0100 << std::endl;
//std::cout << "Check b0010: " << b0010 << std::endl;
//std::cout << "Check b0001: " << b0001 << std::endl;
//std::cout << "Check b0011: " << b0011 << std::endl;
//std::cout << "Check b0101: " << b0101 << std::endl;
//std::cout << "Check b1001: " << b1001 << std::endl;
//std::cout << "Check b0110: " << b0110 << std::endl;
//std::cout << "Check b1010: " << b1010 << std::endl;
//std::cout << "Check b1100: " << b1100 << std::endl;
//std::cout << "Check b0002: " << b0002 << std::endl;
//std::cout << "Check b0020: " << b0020 << std::endl;
//std::cout << "Check b0200: " << b0200 << std::endl;
//std::cout << "Check b2000: " << b2000 << std::endl;
//std::cout << "--------------------------------------------------------------------" << std::endl;

  }


  double p(double T, double muB, double muQ, double muS)
	{
		double x = b0000 + b1000*T + b0100*muB + b0010*muS + b0001*muQ
               + b2000*T*T + b0200*muB*muB + b0020*muS*muS + b0002*muQ*muQ
               + b1100*T*muB + b1010*T*muS + b1001*T*muQ
               + b0110*muB*muS + b0101*muB*muQ + b0011*muS*muQ;
    return x*x;
	}
	
	double s(double T, double muB, double muQ, double muS)
	{
		double x = b0000 + b1000*T + b0100*muB + b0010*muS + b0001*muQ
               + b2000*T*T + b0200*muB*muB + b0020*muS*muS + b0002*muQ*muQ
               + b1100*T*muB + b1010*T*muS + b1001*T*muQ
               + b0110*muB*muS + b0101*muB*muQ + b0011*muS*muQ;
    return 2.0*x*(b1000 + muQ*b1001 + muS*b1010 + muB*b1100 + 2.0*T*b2000);
	}

	double B(double T, double muB, double muQ, double muS)
	{
		double x = b0000 + b1000*T + b0100*muB + b0010*muS + b0001*muQ
               + b2000*T*T + b0200*muB*muB + b0020*muS*muS + b0002*muQ*muQ
               + b1100*T*muB + b1010*T*muS + b1001*T*muQ
               + b0110*muB*muS + b0101*muB*muQ + b0011*muS*muQ;
    return 2.0*x*(b0100 + muQ*b0101 + muS*b0110 + 2.0*muB*b0200 + T*b1100);
	}

	double S(double T, double muB, double muQ, double muS)
	{
		double x = b0000 + b1000*T + b0100*muB + b0010*muS + b0001*muQ
               + b2000*T*T + b0200*muB*muB + b0020*muS*muS + b0002*muQ*muQ
               + b1100*T*muB + b1010*T*muS + b1001*T*muQ
               + b0110*muB*muS + b0101*muB*muQ + b0011*muS*muQ;
    return 2.0*x*(b0010 + muQ*b0011 + 2.0*muS*b0020 + muB*b0110 + T*b1010);
	}

	double Q(double T, double muB, double muQ, double muS)
	{
		double x = b0000 + b1000*T + b0100*muB + b0010*muS + b0001*muQ
               + b2000*T*T + b0200*muB*muB + b0020*muS*muS + b0002*muQ*muQ
               + b1100*T*muB + b1010*T*muS + b1001*T*muQ
               + b0110*muB*muS + b0101*muB*muQ + b0011*muS*muQ;
    return 2.0*x*(b0001 + 2.0*muQ*b0002 + muS*b0011 + muB*b0101 + T*b1001);
	}

  // diagonal susceptibilities
	double P2T2(double T, double muB, double muQ, double muS)
	{
		double x = b0000 + b1000*T + b0100*muB + b0010*muS + b0001*muQ
               + b2000*T*T + b0200*muB*muB + b0020*muS*muS + b0002*muQ*muQ
               + b1100*T*muB + b1010*T*muS + b1001*T*muQ
               + b0110*muB*muS + b0101*muB*muQ + b0011*muS*muQ;
    double yT = b1000 + 2.0*T*b2000 + muB*b1100 + muS*b1010 + muQ*b1001;
		return 2.0*(2.0*b2000*x + yT*yT);
	}
	double P2B2(double T, double muB, double muQ, double muS)
	{
		double x = b0000 + b1000*T + b0100*muB + b0010*muS + b0001*muQ
               + b2000*T*T + b0200*muB*muB + b0020*muS*muS + b0002*muQ*muQ
               + b1100*T*muB + b1010*T*muS + b1001*T*muQ
               + b0110*muB*muS + b0101*muB*muQ + b0011*muS*muQ;
    double yB = b0100 + T*b1100 + 2.0*muB*b0200 + muS*b0110 + muQ*b0101;
		return 2.0*(2.0*b0200*x + yB*yB);
	}
	double P2S2(double T, double muB, double muQ, double muS)
	{
		double x = b0000 + b1000*T + b0100*muB + b0010*muS + b0001*muQ
               + b2000*T*T + b0200*muB*muB + b0020*muS*muS + b0002*muQ*muQ
               + b1100*T*muB + b1010*T*muS + b1001*T*muQ
               + b0110*muB*muS + b0101*muB*muQ + b0011*muS*muQ;
    double yS = b0010 + T*b1010 + muB*b0110 + 2.0*muS*b0020 + muQ*b0011;
		return 2.0*(2.0*b0020*x + yS*yS);
	}
	double P2Q2(double T, double muB, double muQ, double muS)
	{
		double x = b0000 + b1000*T + b0100*muB + b0010*muS + b0001*muQ
               + b2000*T*T + b0200*muB*muB + b0020*muS*muS + b0002*muQ*muQ
               + b1100*T*muB + b1010*T*muS + b1001*T*muQ
               + b0110*muB*muS + b0101*muB*muQ + b0011*muS*muQ;
    double yQ = b0001 + T*b1001 + muB*b0101 + muS*b0011 + 2.0*muQ*b0002;
		return 2.0*(2.0*b0002*x + yQ*yQ);
	}
	
  // off-diagonal susceptibilities
	double P2BQ(double T, double muB, double muQ, double muS)
  {
		double x = b0000 + b1000*T + b0100*muB + b0010*muS + b0001*muQ
               + b2000*T*T + b0200*muB*muB + b0020*muS*muS + b0002*muQ*muQ
               + b1100*T*muB + b1010*T*muS + b1001*T*muQ
               + b0110*muB*muS + b0101*muB*muQ + b0011*muS*muQ;
    double yB = b0100 + T*b1100 + 2.0*muB*b0200 + muS*b0110 + muQ*b0101;
    double yQ = b0001 + T*b1001 + muB*b0101 + muS*b0011 + 2.0*muQ*b0002;
    return 2.0*(b0101*x+yB*yQ);
  }
	double P2BS(double T, double muB, double muQ, double muS)
  {
		double x = b0000 + b1000*T + b0100*muB + b0010*muS + b0001*muQ
               + b2000*T*T + b0200*muB*muB + b0020*muS*muS + b0002*muQ*muQ
               + b1100*T*muB + b1010*T*muS + b1001*T*muQ
               + b0110*muB*muS + b0101*muB*muQ + b0011*muS*muQ;
    double yB = b0100 + T*b1100 + 2.0*muB*b0200 + muS*b0110 + muQ*b0101;
    double yS = b0010 + T*b1010 + muB*b0110 + 2.0*muS*b0020 + muQ*b0011;
    return 2.0*(b0110*x+yB*yS);
  }
	double P2QS(double T, double muB, double muQ, double muS)
  {
		double x = b0000 + b1000*T + b0100*muB + b0010*muS + b0001*muQ
               + b2000*T*T + b0200*muB*muB + b0020*muS*muS + b0002*muQ*muQ
               + b1100*T*muB + b1010*T*muS + b1001*T*muQ
               + b0110*muB*muS + b0101*muB*muQ + b0011*muS*muQ;
    double yS = b0010 + T*b1010 + muB*b0110 + 2.0*muS*b0020 + muQ*b0011;
    double yQ = b0001 + T*b1001 + muB*b0101 + muS*b0011 + 2.0*muQ*b0002;
    return 2.0*(b0011*x+yS*yQ);
  }
	
	double P2TB(double T, double muB, double muQ, double muS)
  {
		double x = b0000 + b1000*T + b0100*muB + b0010*muS + b0001*muQ
               + b2000*T*T + b0200*muB*muB + b0020*muS*muS + b0002*muQ*muQ
               + b1100*T*muB + b1010*T*muS + b1001*T*muQ
               + b0110*muB*muS + b0101*muB*muQ + b0011*muS*muQ;
    double yT = b1000 + 2.0*T*b2000 + muB*b1100 + muS*b1010 + muQ*b1001;
    double yB = b0100 + T*b1100 + 2.0*muB*b0200 + muS*b0110 + muQ*b0101;
    return 2.0*(b1100*x+yT*yB);
  }
	double P2TS(double T, double muB, double muQ, double muS)
  {
		double x = b0000 + b1000*T + b0100*muB + b0010*muS + b0001*muQ
               + b2000*T*T + b0200*muB*muB + b0020*muS*muS + b0002*muQ*muQ
               + b1100*T*muB + b1010*T*muS + b1001*T*muQ
               + b0110*muB*muS + b0101*muB*muQ + b0011*muS*muQ;
    double yT = b1000 + 2.0*T*b2000 + muB*b1100 + muS*b1010 + muQ*b1001;
    double yS = b0010 + T*b1010 + muB*b0110 + 2.0*muS*b0020 + muQ*b0011;
    return 2.0*(b1010*x+yT*yS);
  }
	double P2TQ(double T, double muB, double muQ, double muS)
  {
		double x = b0000 + b1000*T + b0100*muB + b0010*muS + b0001*muQ
               + b2000*T*T + b0200*muB*muB + b0020*muS*muS + b0002*muQ*muQ
               + b1100*T*muB + b1010*T*muS + b1001*T*muQ
               + b0110*muB*muS + b0101*muB*muQ + b0011*muS*muQ;
    double yT = b1000 + 2.0*T*b2000 + muB*b1100 + muS*b1010 + muQ*b1001;
    double yQ = b0001 + T*b1001 + muB*b0101 + muS*b0011 + 2.0*muQ*b0002;
    return 2.0*(b1001*x+yT*yQ);
  }


  void get_full_thermo( const double point[], double results[] )
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

    const double hc = constants::hbarc_MeVfm;
    double T = point[0], muB = point[1], muQ = point[2], muS = point[3];

    double C1T = results[1]/(T*T*T);
    double C1B = results[2]/(T*T*T);
    double C1S = results[3]/(T*T*T);
    double C1Q = results[4]/(T*T*T);

    double C2B2 = results[7]/(T*T);
    double C2Q2 = results[8]/(T*T);
    double C2S2 = results[9]/(T*T);
    double C2BQ = results[10]/(T*T);
    double C2BS = results[11]/(T*T);
    double C2QS = results[12]/(T*T);
    double C2TB = results[13]/(T*T);
    double C2TQ = results[14]/(T*T);
    double C2TS = results[15]/(T*T);
    double C2T2 = results[16]/(T*T);

//    cout << "Check input(C++): " << T << "   " << muB << "   " << muQ << "   " << muS << "   "
//          << C1T << "   " << C1B << "   " << C1S << "   " << C1Q << "   "
//          << C2B2 << "   " << C2Q2 << "   " << C2S2 << "   "
//          << C2BQ << "   " << C2BS << "   " << C2QS << "   "
//          << C2TB << "   " << C2TQ << "   " << C2TS << "   " << C2T2 << endl;

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



  void get_eBSQ( const double point[], double results[] )
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

    results[0]  = eOut;
    results[1]  = BOut;
    results[2]  = SOut;
    results[3]  = QOut;
  }


  void get_sBSQ( const double point[], double results[] )
  {
    // point: (T, muB, muQ, muS)
    const double Tsol   = point[0], muBsol = point[1],
                 muQsol = point[2], muSsol = point[3];

    results[0]  = s(Tsol, muBsol, muQsol, muSsol);
    results[1]  = B(Tsol, muBsol, muQsol, muSsol);
    results[2]  = S(Tsol, muBsol, muQsol, muSsol);
    results[3]  = Q(Tsol, muBsol, muQsol, muSsol);
  }


  //////////////////////////////////////////////////////////////////////////////
  // obtain non-conformal extension
  // results MUST HAVE LENGTH 17
  void get_nonconformal_extension( const double point[],
                                   const double point_projected[],
                                   double results[] )
  {
    // determine parameters from thermodynamic quantities
    set_coeffs( point_projected, results );

    // evaluate extension and return result (stored in results)
    get_full_thermo( point, results );
  }

  //////////////////////////////////////////////////////////////////////////////
  // projects phase diagram point to grid boundary
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