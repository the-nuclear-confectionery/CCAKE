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
    bool print_this_particle         = false;

    int ID                           = -1;
    int btrack                       = 0;
    int Freeze                       = 0;

    double efcheck                   = 0.0;
    double contribution_to_total_E   = 0.0;
    double contribution_to_total_dEz = 0.0;
    double contribution_to_total_Ez  = 0.0;

    Vector<double,2> r;       // transverse position

    //==========================================================================
    // different combinations of densities
    densities input     = {}; // these densities are read in initially and
                              // have physical units ~1/fm^3;
                              // entropy component previously called "s_an"
                              //================================================
    densities smoothed  = {}; // these are the smoothed (propagated) densities
                              // which have units ~/1/fm^2;
                              // entropy component previously called "eta"
                              //================================================
    densities specific  = {}; // these are the densities "per particle" which
                              // are effectively dimensionless;
                              // entropy component previously called "eta_sigma"
                              //================================================
    densities d_dt_spec = {}; // these are the TIME DERIVATIVES of the specific
                              // densities above;
                              // entropy component previously called "detasigma_dt"
                              //================================================
    densities norm_spec = {}; // gives the normalizations of each of the
                              // specific densities above, can choose different
                              // values for different densities by convenience;
                              // entropy component previously called "sigmaweight"
                              //================================================

    // structs for hydrodynamic and thermodynamic information
    hydrodynamic_info  hydro  = {};
    thermodynamic_info thermo = {};


    //==========================================================================
    // FUNCTIONS AND ROUTINES

    // Constructors and destructors.
    Particle();
    
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
