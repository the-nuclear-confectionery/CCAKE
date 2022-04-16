#ifndef PARTICLE_H
#define PARTICLE_H

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "eos.h"
#include "matrix.h"
#include "vector.h"

class Particle
{
  public:

    // Constructors and destructors.
    Particle();
    Particle(vector<double> &fields);
    
    // copy-constructor
    Particle( const Particle& p );
   ~Particle(){}

  bool operator==( const Particle & ) const;

  // use this to set equation of state object before creating particles
  void set_EquationOfStatePtr( EquationOfState * eosPtr_in );
  void set_SettingsPtr(Settings * settingsPtr_in);

  double locate_phase_diagram_point_eBSQ( double e_In );
  double locate_phase_diagram_point_eBSQ( double e_In, double rhoB_In,
                                          double rhoS_In, double rhoQ_In );
  void locate_phase_diagram_point_sBSQ( double s_In );
  void locate_phase_diagram_point_sBSQ( double s_In, double rhoB_In,
                                        double rhoS_In, double rhoQ_In );

  //private:

    // thermodynamic quantities (with default initialization)
    struct thermodynamic_info
    {
      // (T,mu_i) coordinates depend on which EoS was used!
      string eos_name = "";

      double T    = 0.0, muB  = 0.0, muS  = 0.0, muQ  = 0.0;
      double e    = 0.0, s    = 0.0, rhoB = 0.0, rhoS = 0.0, rhoQ = 0.0,
             p    = 0.0, cs2  = 0.0, w    = 0.0, A    = 0.0;
      double dwds = 0.0, dwdB = 0.0, dwdS = 0.0, dwdQ = 0.0;

      public:
        void set(EquationOfState & eos)
        {
          eos_name = eos.get_current_eos_name();

          T    = eos.T();
          muB  = eos.muB();
          muS  = eos.muS();
          muQ  = eos.muQ();

          p    = eos.p();
          s    = eos.s();
          rhoB = eos.B();
          rhoS = eos.S();
          rhoQ = eos.Q();
          e    = eos.e();
          w    = eos.w();
          A    = eos.A();
          cs2  = eos.cs2();
          dwds = eos.dwds();
          dwdB = eos.dwdB();
          dwdS = eos.dwdS();
          dwdQ = eos.dwdQ();
        }
    };
    thermodynamic_info thermo = {};


//  public:

    EquationOfState * eosPtr = nullptr;
    Settings * settingsPtr   = nullptr;

    // getter functions
    double T()    { return thermo.T;    }
    double muB()  { return thermo.muB;  }
    double muS()  { return thermo.muS;  }
    double muQ()  { return thermo.muQ;  }

    double p()    { return thermo.p;    }
    double s()    { return thermo.s;    }
    double e()    { return thermo.e;    }
    double rhoB() { return thermo.rhoB; }
    double rhoS() { return thermo.rhoS; }
    double rhoQ() { return thermo.rhoQ; }
    double w()    { return thermo.w;    }
    double A()    { return thermo.A;    }
    double cs2()  { return thermo.cs2;  }
    double dwds() { return thermo.dwds; }
    double dwdB() { return thermo.dwdB; }
    double dwdS() { return thermo.dwdS; }
    double dwdQ() { return thermo.dwdQ; }

    string get_current_eos_name() { return thermo.eos_name; }


  // rename these functions and their arguments
  void frzcheck( double tin, int &count, int N );
  void calc(double tin);
  void calcbsq(double tin);
  void return_bsqsv_A();
  void bsqsvsigset(double tin, int i);
  void setvisc( int etaconst, double bvf, double svf, double zTc, double sTc, double sig, int type );
  void sets(double tin2, bool is_first_timestep);
  void setvar();
  double gamcalc();
  double Bsub();

  // members
  int btrack             = 0;
  //int count              = 0;
  int Freeze             = 0;

  double Agam            = 0.0;
  double Agam2           = 0.0;
  double sigmaweight     = 0.0; // specific volume per particle (times s_an)
  double rhoB_weight     = 0.0; // specific volume per particle (without rhoB_an)
  double rhoS_weight     = 0.0; // specific volume per particle (without rhoS_an)
  double rhoQ_weight     = 0.0; // specific volume per particle (without rhoQ_an)
  double transverse_area = 0.0; // dx * dy
  double ets1            = 0.0;
  double ets2            = 0.0;
  double ets3            = 0.0;
  double ets4            = 0.0;
  double b1              = 0.0;
  double b2              = 0.0;
  double b3              = 0.0;
  double b4              = 0.0;
  double bn1             = 0.0;
  double bn2             = 0.0;
  double bn3             = 0.0;
  double bn4             = 0.0;
  double shv33           = 0.0;

  double div_u           = 0.0; // four-divergence of relativistic velocity
  double gamma           = 0.0; // Lorentz factor
  double s_sub           = 0.0;
  double e_sub           = 0.0;
  double s_an            = 0.0;
  double s_rat           = 0.0;
  double sigsub          = 0.0;
  double eta_sigma       = 0.0; // Ratio entropy/especific volume
  double detasigma_dt    = 0.0;
  double Bulk            = 0.0; // Bulk Viscosity weight
  double bigPI           = 0.0; // total bulk viscosity
  double C               = 0.0;
  double tauRelax        = 0.0; // Bulk Relaxation time
  double stauRelax       = 0.0; // Shear Relxation time
  double dBulk_dt        = 0.0; // derivative Bulk Viscosity
  double zeta            = 0.0; // bulk coefficient
  double setas           = 0.0; 
  double Ctot            = 0.0;
  double Btot            = 0.0;

  double sv_eta          = 0.0;
  double taupi           = 0.0;

  ////////////////////////////////////////////////////////////////////////////
  //                         Fluid Variables                                //
  ////////////////////////////////////////////////////////////////////////////
  double sigma           = 0.0; // especific volume
  double dsigma_dt       = 0.0; // derivative of especific volume

  double dw_ds           = 0.0; // derivative of the enthalpy on entropy
  double eta             = 0.0; // entropy density
  double rhoB_an         = 0.0;
  double rhoB_sub        = 0.0; // Baryon density
  double rhoS_an         = 0.0;
  double rhoS_sub        = 0.0; // strange density
  double rhoQ_an         = 0.0;
  double rhoQ_sub        = 0.0; // electric charge density
  double eden            = 0.0;

  double bigtheta        = 0.0;

  double B               = 0.0; // Baryon density
  double S               = 0.0; // Baryon density
  double Q               = 0.0; // baryon, strange, electric charge
  double drhoB_dt        = 0.0; // Baryon density
  double drhoS_dt        = 0.0; // Baryon density
  double drhoQ_dt        = 0.0;

  double g2              = 0.0;
  double g3              = 0.0;
  double gt              = 0.0;
  double eta_o_tau       = 0.0;
  double dwdsT1          = 0.0;
  double sigl            = 0.0;
  double bcheck          = 0.0;
  double check           = 0.0;
  double saves           = 0.0;
  double inside          = 0.0;

  double freezeoutT      = 0.0;
  double efcheck         = 0.0;

  // vector members
  Vector<double,2> r;                                       // position
  Vector<double,2> v;                                       // velocity
  Vector<double,2> u;                                       // relativistic velocity
  Vector<double,2> qmom;
  Vector<double,2> du_dt, gradsig;                          // relativistic velocity derivative
  Vector<double,2> k1, k2, k3, k4, ksub;
  Vector<double,2> r1, r2, r3, r4, rsub;

  Vector<double,2> gradP;                                   // Gradient of Pressure
  Vector<double,2> gradBulk, gradrhoB, gradrhoS, gradrhoQ;  // Gradient of Bulk Viscosity
  Vector<double,2> gradsigma;                               // Gradient of especific volume
  Vector<double,2> divshear, gradshear;


  // matrix members
  Matrix<double,2,2> Msub(int i);
  Matrix<double,2,2> dpidtsub();
  Matrix<double,2,2> Imat;
  Matrix<double,2,2> gradV, gradU;                          // Gradient of velocity needed for shear
  Matrix<double,2,2> /*shv0, */dshv_dt;
  Matrix<double,2,2> shv1, shv2, shv3, svh4;
  Matrix<double,2,2> uu, pimin, piu, piutot;
  Matrix<double,3,3> shv;


  // freeze out struct thingy?
  struct FRZ
  {
      double t = 0.0, s = 0.0, e = 0.0, rhoB = 0.0, rhoS = 0.0, rhoQ = 0.0,
             T = 0.0, muB = 0.0, muS = 0.0, muQ = 0.0, theta = 0.0, bulk = 0.0,
             sigma = 0.0, shear33 = 0.0, inside = 0.0;
      Vector<double,2> r, u, gradP;
      Matrix<double,3,3> shear;
  };

  FRZ frz1   = {};
  FRZ frz2   = {};
  FRZ fback  = {};
  FRZ fback2 = {};
  FRZ fback3 = {};
  FRZ fback4 = {};

////////////////////////////////////////////////////////////////////////////////
// NEW VARIABLES AND FUNCTIONS TO CLEAN UP EQUATIONS_OF_MOTION BELOW THIS LINE
////////////////////////////////////////////////////////////////////////////////
double vartheta = 0.0;

void set_vartheta(double t);
double get_Theta_force();
double get_Pi_force(double t);
double get_aleph_force(double t);
Vector<double,2> get_Theta_mass();
Vector<double,2> get_Pi_mass();
Vector<double,2> get_aleph_mass();






};

#endif
