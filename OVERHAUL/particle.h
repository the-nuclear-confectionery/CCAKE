#ifndef PARTICLE_H
#define PARTICLE_H

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
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
    //Particle( const Particle& p );
   ~Particle(){}

  bool operator==( const Particle & ) const;

  // use this to set equation of state object before creating particles
  void set_EquationOfStatePtr( EquationOfState * eosPtr_in );

  double locate_phase_diagram_point_eBSQ( double e_In );
  double locate_phase_diagram_point_eBSQ( double e_In, double rhoB_In,
                                          double rhoS_In, double rhoQ_In );
  void locate_phase_diagram_point_sBSQ( double s_In );
  void locate_phase_diagram_point_sBSQ( double s_In, double rhoB_In,
                                        double rhoS_In, double rhoQ_In );

  //private:
    // kinematic quantities
//    vector<double> r;  // (x,y) in fm
//    vector<double> v;  // velocity
//    vector<double> u;  // relativistic velocity

    // thermodynamic quantities (with default initialization)
    struct thermodynamic_info
    {
      double T    = 0.0, muB  = 0.0, muS  = 0.0, muQ  = 0.0;
      double e    = 0.0, s    = 0.0, rhoB = 0.0, rhoS = 0.0, rhoQ = 0.0,
             p    = 0.0, cs2  = 0.0, w    = 0.0, A    = 0.0;
      double dwds = 0.0, dwdB = 0.0, dwdS = 0.0, dwdQ = 0.0;

      public:
        void set(EquationOfState & eos)
        {
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
    thermodynamic_info thermo;


//  public:

    //static EquationOfState eos;	//use one copy of EOS for all particles
    EquationOfState * eosPtr;

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


    // rename these functions and their arguments
//    void calcbsq( double tin );
//    void bsqsvsigset( double tin, int i );
//    void return_bsqsv_A();
//    void setvisc( int etaconst, double bvf, double svf, double zTc, double sTc,
//                  double sig, int type );
//    void frzcheck( double tin, int & count, int N );
  double gamcalc();
  void frzcheck( double tin, int &count, int N );
  void calc(double tin);
  void calcbsq(double tin);
  void return_bsqsv_A();
  double Bsub();
  Matrix<double,2,2> Msub(int i);
  Matrix<double,2,2> dpidtsub();
  void bsqsvsigset(double tin,int i);
  void setvisc( int etaconst, double bvf, double svf, double zTc,
                          double sTc, double sig, int type );
  void sets(double tin2);
  void setvar();















  // old Particle class variables

    Matrix <double,2,2> Imat;


    int btrack;
    double Agam, Agam2;
    double sigmaweight;        // specific volume per particle (times s_an)
    double rho_weight;		   // specific volume per particle (without s_an)
    double transverse_area;	   // dx * dy
    Vector<double,2> r;                   // position
    Vector<double,2> v;                   // velocity
    Vector<double,2> u;                   // relativistic velocity
    Vector<double,2> qmom;
    Vector<double,2> du_dt,gradsig;               // relativistic velocity derivative
    Vector<double,2> k1,k2,k3,k4,ksub;
    Vector<double,2> r1,r2,r3,r4,rsub;
    double ets1,ets2,ets3,ets4;
    double b1,b2,b3,b4;
    double bn1,bn2,bn3,bn4;
    Matrix<double,2,2> shv1,shv2,shv3,svh4;
    double shv33;


    double freezeoutT;

    struct FRZ
    {
        Vector<double,2> r,u,gradP;
        double t,s,T,theta,bulk,sigma,shear33,inside;
        Matrix<double,3,3> shear;
    };

    FRZ frz1,frz2,fback,fback2,fback3,fback4;

    double div_u;              // four-divergence of relativistic velocity
    double gamma;               // Lorentz factor
    double s_sub,e_sub,s_an,s_rat,sigsub;
    double eta_sigma;          // Ratio entropy/especific volume
    double detasigma_dt;
    double Bulk;               // Bulk Viscosity weight
    double bigPI;        // total bulk viscosity
    double C;
    double tauRelax;           // Bulk Relaxation time
    double stauRelax;        // Shear Relxation time;
    double dBulk_dt;           // derivative Bulk Viscosity
    double zeta;               // bulk coefficient
    double setas;
    int count;
    int Freeze;
    double Ctot, Btot;
    Matrix<double,3,3> shv;
    Matrix<double,2,2> shv0,dshv_dt;
    Vector<double,2>  divshear,gradshear;

    double sv_eta,taupi;

////////////////////////////////////////////////////////////////////////////////
//                           Fluid Variables                                  //
////////////////////////////////////////////////////////////////////////////////
    double sigma;              // especific volume
    double dsigma_dt;          // derivative of especific volume
    Vector<double,2> gradP;              // Gradient of Pressure
    Matrix<double,2,2> gradV,gradU;          // Gradient of velocity needed for shear

    Vector<double,2> gradBulk,gradrhoB,gradrhoS,gradrhoQ;           // Gradient of Bulk Viscosity
    Vector<double,2> gradsigma;          // Gradient of especific volume

    double dw_ds;              // derivative of the enthalpy on entropy
    double eta;                   // entropy density
    double rhoB_an, rhoB_sub;                   // Baryon density
    double rhoS_an, rhoS_sub;                   // strange density
    double rhoQ_an, rhoQ_sub;                   // electric charge density
    double eden;

    double bigtheta;

    double B, S, Q;						// baryon, strange, electric charge
    double drhoB_dt,drhoS_dt,drhoQ_dt;

    double g2,g3,gt;
    double eta_o_tau,dwdsT1,sigl;
    Matrix <double,2,2> uu,pimin,piu,piutot;
    double bcheck;
    double check;
    double saves;
    double inside;





};

#endif
