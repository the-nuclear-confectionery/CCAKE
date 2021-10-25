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

class Particle
{

  public:

    // Constructors and destructors.
    Particle();
   ~Particle();

  // use this to set equation of state object before creating particles
  static void set_equation_of_state( EquationOfState & eos );

  double locate_phase_diagram_point_eBSQ( double e_In );
  double locate_phase_diagram_point_eBSQ( double e_In, double rhoB_In,
                                          double rhoS_In, double rhoQ_In );
  void locate_phase_diagram_point_sBSQ( double s_In );
  void locate_phase_diagram_point_sBSQ( double s_In, double rhoB_In,
                                        double rhoS_In, double rhoQ_In );

  private:
    // kinematic quantities
    vector<double> r;  // (x,y) in fm
    vector<double> v;  // velocity
    vector<double> u;  // relativistic velocity

    // thermodynamic quantities (with default initialization)
    struct thermodynamic_info
    {
      double T    = 0.0, muB  = 0.0, muS  = 0.0, muQ  = 0.0;
      double e    = 0.0, s    = 0.0, B    = 0.0, S    = 0.0, Q = 0.0,
             p    = 0.0, cs2  = 0.0, w    = 0.0, A    = 0.0;
      double dwds = 0.0, dwdB = 0.0, dwdS = 0.0, dwdQ = 0.0;

      public:
        void set(const EquationOfState & eos)
        {
          T    = eos.T();
          muB  = eos.muB();
          muS  = eos.muS();
          muQ  = eos.muQ();

          p    = eos.p();
          s    = eos.s();
          B    = eos.B();
          S    = eos.S();
          Q    = eos.Q();
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

    static EquationOfState eos;	//use one copy of EOS for all particles
    //EquationOfState * eosPtr = nullptr;

  public:

    // getter functions
    double T()    { return thermo.T;    }
    double muB()  { return thermo.muB;  }
    double muS()  { return thermo.muS;  }
    double muQ()  { return thermo.muQ;  }

    double p()    { return thermo.p;    }
    double s()    { return thermo.s;    }
    double e()    { return thermo.e;    }
    double B()    { return thermo.B;    }
    double S()    { return thermo.S;    }
    double Q()    { return thermo.Q;    }
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


};

#endif