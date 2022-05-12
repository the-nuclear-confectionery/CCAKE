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

#include "constants.h"
#include "densities.h"
#include "hydrodynamic_info.h"
#include "matrix.h"
#include "settings.h"
#include "thermodynamic_info.h"
#include "vector.h"

class Particle
{
  private:

    static constexpr int VERBOSE = 0;

    Settings * settingsPtr   = nullptr;

  public:

    //==========================================================================
    // MEMBERS
    bool print_this_particle = false;

    int ID                 = -1;
    int btrack             = 0;
    int Freeze             = 0;

    double contribution_to_total_E   = 0.0;
    double contribution_to_total_dEz = 0.0;
    double contribution_to_total_Ez  = 0.0;

    double ets1 = 0.0, ets2 = 0.0, ets3 = 0.0, ets4 = 0.0;
    double b1 = 0.0, b2 = 0.0, b3 = 0.0, b4 = 0.0;
    Vector<double,2> k1, k2, k3, k4;
    Vector<double,2> r1, r2, r3, r4;

    Vector<double,2> r;       // transverse position

    // different combinations of densities
    densities input    = {};  // these densities are read in initially which
                              // have physical units ~1/fm^3
                              // entropy component previously called "s_an"
    densities smoothed = {};  // these are the smoothed (propagated) densities
                              // which have units ~/1/fm^2
                              // entropy component previously called "eta"
    densities specific = {};  // these are the densities "per particle" which
                              // are effectively dimensionless
                              // entropy component previously called "eta_sigma"


    ////////////////////////////////////////////////////////////////////////////
    //                         Fluid Variables                                //
    ////////////////////////////////////////////////////////////////////////////
    //double e_sub           = 0.0;
    //double s_an            = 0.0;
    //double eta             = 0.0; // entropy density
    //double rhoB_an         = 0.0;
    //double rhoB_sub        = 0.0; // Baryon density
    //double rhoS_an         = 0.0;
    //double rhoS_sub        = 0.0; // strange density
    //double rhoQ_an         = 0.0;
    //double rhoQ_sub        = 0.0; // electric charge density

    double eta_sigma       = 0.0; // Ratio entropy/especific volume
    double B               = 0.0; // Baryon per particle
    double S               = 0.0; // Net strangeness per particle
    double Q               = 0.0; // Electric charge per particle

    double sigmaweight     = 0.0; // specific volume per particle (times s_an)
    double rhoB_weight     = 0.0; // specific volume per particle (without rhoB_an)
    double rhoS_weight     = 0.0; // specific volume per particle (without rhoS_an)
    double rhoQ_weight     = 0.0; // specific volume per particle (without rhoQ_an)


    hydrodynamic_info  hydro  = {};
    thermodynamic_info thermo = {};





    // freeze out struct
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

    double efcheck = 0.0;





    //==========================================================================
    // FUNCTIONS AND ROUTINES

    // Constructors and destructors.
    Particle();
    Particle(vector<double> &fields);
    
    // copy-constructor
    Particle( const Particle& p );
   ~Particle(){}

    bool operator==( const Particle & ) const;

    // use this to set equation of state object before creating particles
    void set_SettingsPtr(Settings * settingsPtr_in) { settingsPtr = settingsPtr_in; }

    void evaluate_time_derivatives( double t );

    //==========================================================================
    // getter functions for thermodynamic information
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

    //==========================================================================
    // getter functions for hydrodynamic information
    //...
    //...add these later...
    //...



    // rename these functions and their arguments
    void reset_pi_tensor(double tin2);
    double gamcalc();

};

#endif
