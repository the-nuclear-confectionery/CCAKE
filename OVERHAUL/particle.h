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

class particle
{
  // Constructors and destructors.
  // particle();
  // ~particle();




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
        void set(const eos & EOS)
        {
          T    = EOS.T();
          muB  = EOS.muB();
          muS  = EOS.muS();
          muQ  = EOS.muQ();

          p    = EOS.p();
          s    = EOS.s();
          B    = EOS.B();
          S    = EOS.S();
          Q    = EOS.Q();
          e    = EOS.e();
          w    = EOS.w();
          A    = EOS.A();
          cs2  = EOS.cs2();
          dwds = EOS.dwds();
          dwdB = EOS.dwdB();
          dwdS = EOS.dwdS();
          dwdQ = EOS.dwdQ();
        }
    }
    thermodynamic_info thermo;

    static eos EOS;	//use one copy of EOS for all particles


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


}

#endif