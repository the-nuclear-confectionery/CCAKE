#ifndef SYSTEM_STATE_H
#define SYSTEM_STATE_H

#include <iostream>
#include <string>
#include <vector>

#include "eos.h"
#include "kernel.h"
#include "linklist.h"
#include "matrix.h"
#include "particle.h"
#include "settings.h"
#include "vector.h"

using std::string;
using std::vector;

class SystemState
{
  friend class InputOutput;
  friend class SPHWorkstation;
  friend class Evolver;
  friend class FreezeOut;

  public:

    SystemState(){}
    ~SystemState(){}

    void set_EquationOfStatePtr( EquationOfState * eosPtr_in );
    void set_SettingsPtr( Settings * settingsPtr_in );

    bool do_freeze_out = true;

    double h          = 0.0;
    double t          = 0.0;
    double dt         = 0.0;

    double S          = 0.0;
    double S0         = 0.0;
    double Btotal     = 0.0;
    double Btotal0    = 0.0;
    double Stotal     = 0.0;
    double Stotal0    = 0.0;
    double Qtotal     = 0.0;
    double Qtotal0    = 0.0;

    double E          = 0.0;
    double Ez         = 0.0;
    double E0         = 0.0;
    double Eloss      = 0.0;
    double dEz        = 0.0;
    double Etot       = 0.0;

    double efcheck    = 0.0;
    double sfcheck    = 0.0;

   /* INTS */
    int number_part = 0;
    int n_particles = 0;
    int rk2         = 0;

    int number_of_elapsed_timesteps = 0;


  private:
    EquationOfState * eosPtr        = nullptr;
    Settings        * settingsPtr   = nullptr;


    // the vector of particles
    vector< Particle > particles;

    // the linklist specifying nearest neighbors
    LinkList linklist; 

    // used to track when a particle fails to find a solution in any EoS
    vector<int> particles_out_of_grid;

    vector<int> list;

    // eccentricities
    vector<double> timesteps;
    vector<double> e_2_X, e_2_P;

  public:
    // initialize system state, linklist, etc.
    void initialize();
    void initialize_linklist();

    // check conserved quantities
    void conservation_entropy();
    void conservation_energy();
    void conservation_BSQ();

    void compute_eccentricities();
    void compute_e_2_P();
    void compute_e_2_X();

    int n(){ return n_particles; }
    double get_particle_T(int id) {return particles[id].T();}
    double get_particle_Freeze(int id) {return particles[id].Freeze;}
    int get_frozen_out_count()
    {
      int total_frz_out = 0;
      for ( auto & p : particles)
        if (p.Freeze == 4) total_frz_out++;
      return total_frz_out;
    }

};

#endif