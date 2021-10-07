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

    // thermodynamic quantities
    struct thermodynamic_info
    {
      double T, muB, muS, muQ;
      double e, p, s, B, S, Q, cs2, w, A;
      double dwds, dwdB, dwdS, dwdQ;
    }
    thermodynamic_info thermo;

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
    double dwds() { return thermo.dwds; }
    double dwdB() { return thermo.dwdB; }
    double dwdS() { return thermo.dwdS; }
    double dwdQ() { return thermo.dwdQ; }


  public:



}

#endif