#ifndef EOS_TANH_CONFORMAL_H
#define EOS_TANH_CONFORMAL_H

#include <cmath>
#include <vector>

#include "constants.h"
#include "eos_base.h"
#include "eos_header.h"

class EoS_tanh_conformal: public EoS_base
{
private:
  double c, T0, muB0, muS0, muQ0;
  double Ts, Tc;

  double t_factor, tp_factor, tpp_factor;

  inline double pow2( double x ) { return x*x; }
  inline double sech( double x ) { return 1.0/cosh(x); }

  void set_hyp( double T )
  {
    double t_arg = std::tanh( (T-Tc)/Ts );
    t_factor   = 0.5*(1.0+t_arg);
    tp_factor  = 0.5*pow2( sech( (T-Tc)/Ts ) )/Ts;
    // tp_factor  = 0.5*( [](int x){ return x*x; }(sech( (T-Tc)/Ts) )/Ts;
    tpp_factor = -2.0*t_arg*tp_factor/Ts;
  }

public:
  // default constructor/destructor
  EoS_tanh_conformal(){}
  virtual ~EoS_tanh_conformal(){}

  EoS_tanh_conformal( const double c_in, const double T0_in, const double muB0_in,
                 const double muS0_in, const double muQ0_in,
                 const double Tc_in, const double Ts_in,
                 const std::vector<double> & tbqs_minima_in,
                 const std::vector<double> & tbqs_maxima_in,
                 const std::string & name_in = "tanh_conformal")
    { c = c_in; T0 = T0_in; muB0 = muB0_in; muS0 = muS0_in; muQ0 = muQ0_in;
      Tc = Tc_in; Ts = Ts_in;
      tbqs_minima = tbqs_minima_in; tbqs_maxima = tbqs_maxima_in; 
      tbqs_minima_no_ext = tbqs_minima_in; tbqs_maxima_no_ext = tbqs_maxima_in;
      name = name_in; }


  double p(double T, double muB, double muQ, double muS)
	{
		double x = (T/T0)*(T/T0) + (muB/muB0)*(muB/muB0)
                + (muS/muS0)*(muS/muS0) + (muQ/muQ0)*(muQ/muQ0);
    double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
		return cp*x*x;
	}
	
	double s(double T, double muB, double muQ, double muS)
	{
		double x = (T/T0)*(T/T0) + (muB/muB0)*(muB/muB0)
                + (muS/muS0)*(muS/muS0) + (muQ/muQ0)*(muQ/muQ0);
    double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
		return 4.0*cp*T*x/(T0*T0);
	}

	double B(double T, double muB, double muQ, double muS)
	{
		double x = (T/T0)*(T/T0) + (muB/muB0)*(muB/muB0)
                + (muS/muS0)*(muS/muS0) + (muQ/muQ0)*(muQ/muQ0);
    double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
		return 4.0*cp*muB*x/(muB0*muB0);
	}

	double S(double T, double muB, double muQ, double muS)
	{
		double x = (T/T0)*(T/T0) + (muB/muB0)*(muB/muB0)
                + (muS/muS0)*(muS/muS0) + (muQ/muQ0)*(muQ/muQ0);
    double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
		return 4.0*cp*muS*x/(muS0*muS0);
	}

	double Q(double T, double muB, double muQ, double muS)
	{
		double x = (T/T0)*(T/T0) + (muB/muB0)*(muB/muB0)
                + (muS/muS0)*(muS/muS0) + (muQ/muQ0)*(muQ/muQ0);
    double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
		return 4.0*cp*muQ*x/(muQ0*muQ0);
	}

	double P2B2(double T, double muB, double muQ, double muS)
	{
		double x = (T/T0)*(T/T0) + (muB/muB0)*(muB/muB0)
                + (muS/muS0)*(muS/muS0) + (muQ/muQ0)*(muQ/muQ0);
    double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
		return cp*(8.0*muB*muB/(muB0*muB0*muB0*muB0) + 4.0*x/(muB0*muB0));
	}
	double P2Q2(double T, double muB, double muQ, double muS)
	{
		double x = (T/T0)*(T/T0) + (muB/muB0)*(muB/muB0)
                + (muS/muS0)*(muS/muS0) + (muQ/muQ0)*(muQ/muQ0);
    double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
		return cp*(8.0*muQ*muQ/(muQ0*muQ0*muQ0*muQ0) + 4.0*x/(muQ0*muQ0));
	}
	double P2S2(double T, double muB, double muQ, double muS)
	{
		double x = (T/T0)*(T/T0) + (muB/muB0)*(muB/muB0)
                + (muS/muS0)*(muS/muS0) + (muQ/muQ0)*(muQ/muQ0);
    double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
		return cp*(8.0*muS*muS/(muS0*muS0*muS0*muS0) + 4.0*x/(muS0*muS0));
	}
	
	double P2BQ(double T, double muB, double muQ, double muS)
      { double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
        return (8.0*cp*muB*muQ/(muB0*muB0*muQ0*muQ0)); }
	double P2BS(double T, double muB, double muQ, double muS)
			{ double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
        return (8.0*cp*muB*muS/(muB0*muB0*muS0*muS0)); }
	double P2QS(double T, double muB, double muQ, double muS)
			{ double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
        return (8.0*cp*muS*muQ/(muS0*muS0*muQ0*muQ0)); }
	
	double P2TB(double T, double muB, double muQ, double muS)
			{ double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
        return (8.0*cp*T*muB/(T0*T0*muB0*muB0)); }
	double P2TQ(double T, double muB, double muQ, double muS)
			{ double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
        return (8.0*cp*T*muQ/(T0*T0*muQ0*muQ0)); }
	double P2TS(double T, double muB, double muQ, double muS)
			{ double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
        return (8.0*cp*T*muS/(T0*T0*muS0*muS0)); }
	double P2T2(double T, double muB, double muQ, double muS)
	{
    double x = (T/T0)*(T/T0) + (muB/muB0)*(muB/muB0)
                + (muS/muS0)*(muS/muS0) + (muQ/muQ0)*(muQ/muQ0);
    double cp = c*T0*T0*T0*T0;  // c dimensionless, cp - 1/fm^4
		return cp*(8.0*T*T/(T0*T0*T0*T0) + 4.0*x/(T0*T0));
	}

  void get_eBSQ( double point[], double results[] )
  {
    // point: (T, muB, muQ, muS)
    const double Tsol   = point[0], muBsol = point[1],
                 muQsol = point[2], muSsol = point[3];

//cout << "EVALUATE eBSQ at "
//      << Tsol << "   " << muBsol << "   "
//      << muQsol << "   " << muSsol << endl;

    set_hyp( Tsol );

    double POut_c      = p(Tsol, muBsol, muQsol, muSsol);
    double sOut_c      = s(Tsol, muBsol, muQsol, muSsol);
    double BOut_c      = B(Tsol, muBsol, muQsol, muSsol);
    double SOut_c      = S(Tsol, muBsol, muQsol, muSsol);
    double QOut_c      = Q(Tsol, muBsol, muQsol, muSsol);
    double POut        = t_factor*POut_c;
    double sOut        = t_factor*sOut_c + tp_factor*POut_c;
    double BOut        = t_factor*BOut_c;
    double SOut        = t_factor*SOut_c;
    double QOut        = t_factor*QOut_c;
    double eOut        = (sOut*Tsol + muBsol*BOut + muQsol*QOut + muSsol*SOut) - POut;

//cout << "eOut = " << eOut << endl;
//cout << "Misc: " << POut_c << "   " << sOut_c << "   " << BOut_c << "   "
//      << SOut_c << "   " << QOut_c << endl;

    results[0]         = eOut;
    results[1]         = BOut;
    results[2]         = SOut;
    results[3]         = QOut;
  }

  void get_sBSQ( double point[], double results[] )
  {
    // point: (T, muB, muQ, muS)
    const double Tsol   = point[0], muBsol = point[1],
                 muQsol = point[2], muSsol = point[3];

    set_hyp( Tsol );

    double POut_c      = p(Tsol, muBsol, muQsol, muSsol);
    double sOut_c      = s(Tsol, muBsol, muQsol, muSsol);
    double BOut_c      = B(Tsol, muBsol, muQsol, muSsol);
    double SOut_c      = S(Tsol, muBsol, muQsol, muSsol);
    double QOut_c      = Q(Tsol, muBsol, muQsol, muSsol);
    double sOut        = t_factor*sOut_c + tp_factor*POut_c;
    double BOut        = t_factor*BOut_c;
    double SOut        = t_factor*SOut_c;
    double QOut        = t_factor*QOut_c;

    // point: (T, muB, muQ, muS)
    results[0] = t_factor*sOut_c + tp_factor*POut_c;
    results[1] = t_factor*BOut_c;
    results[2] = t_factor*SOut_c;
    results[3] = t_factor*QOut_c;
  }

  void get_full_thermo( double point[], double results[] )
  {
    // point: (T, muB, muQ, muS)
    const double Tsol   = point[0], muBsol = point[1],
                 muQsol = point[2], muSsol = point[3];

    set_hyp( Tsol );

//cout << "EVALUATE full_thermo at "
//      << Tsol << "   " << muBsol << "   "
//      << muQsol << "   " << muSsol << endl;

    double POut_c      = p(Tsol, muBsol, muQsol, muSsol);
    double sOut_c      = s(Tsol, muBsol, muQsol, muSsol);
    double BOut_c      = B(Tsol, muBsol, muQsol, muSsol);
    double SOut_c      = S(Tsol, muBsol, muQsol, muSsol);
    double QOut_c      = Q(Tsol, muBsol, muQsol, muSsol);
    double POut        = t_factor*POut_c;
    double sOut        = t_factor*sOut_c + tp_factor*POut_c;
    double BOut        = t_factor*BOut_c;
    double SOut        = t_factor*SOut_c;
    double QOut        = t_factor*QOut_c;
    double eOut = (sOut*Tsol + muBsol*BOut + muQsol*QOut + muSsol*SOut) - POut;

//cout << "eOut = " << eOut << endl;
//cout << "Misc: " << POut_c << "   " << sOut_c << "   " << BOut_c << "   "
//      << SOut_c << "   " << QOut_c << endl;

    results[0]  = POut;
    results[1]  = sOut;
    results[2]  = BOut;
    results[3]  = SOut;
    results[4]  = QOut;
    results[5]  = eOut;
    //results[6]  = 1.0/3.0;  // conformal

    results[7]  = t_factor*P2B2(Tsol, muBsol, muQsol, muSsol);
    results[8]  = t_factor*P2Q2(Tsol, muBsol, muQsol, muSsol);
    results[9]  = t_factor*P2S2(Tsol, muBsol, muQsol, muSsol);
    results[10] = t_factor*P2BQ(Tsol, muBsol, muQsol, muSsol);
    results[11] = t_factor*P2BS(Tsol, muBsol, muQsol, muSsol);
    results[12] = t_factor*P2QS(Tsol, muBsol, muQsol, muSsol);
    results[13] = t_factor*P2TB(Tsol, muBsol, muQsol, muSsol) + tp_factor*BOut_c;
    results[14] = t_factor*P2TQ(Tsol, muBsol, muQsol, muSsol) + tp_factor*QOut_c;
    results[15] = t_factor*P2TS(Tsol, muBsol, muQsol, muSsol) + tp_factor*SOut_c;
    results[16] = t_factor*P2T2(Tsol, muBsol, muQsol, muSsol)
                   + 2.0*tp_factor*sOut_c + tpp_factor*POut_c;

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
  
};

#endif