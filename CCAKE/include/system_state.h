#ifndef SYSTEM_STATE_H
#define SYSTEM_STATE_H

#include <Cabana_Core.hpp>

#include <iostream>
#include <string>
#include <vector>

#include "constants.h"
#include "eos.h"
#include "formatted_output.h"
#include "kernel.h"
#include "matrix.h"
#include "particle.h"
#include "settings.h"
#include "stopwatch.h"
#include "vector.h"
#include "utilities.h"
#include "densities.h"

using std::string;
using std::vector;


namespace ccake{
template <unsigned int D>
class SystemState
{

  public:

    SystemState() = delete; ///< Default constructor is deleted. Settings must be passed in.
    SystemState( shared_ptr<Settings> settingsPtr_in): settingsPtr(settingsPtr_in) {};
     ~SystemState(){}

    void set_EquationOfStatePtr( EquationOfState * eosPtr_in );
    void set_SettingsPtr( Settings * settingsPtr_in );

    bool do_freeze_out = true;

    double hT          = 0.0;
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


    std::vector<Particle<D>> particles;     ///< Vector of particles
    Cabana::AoSoA<CabanaParticle, DeviceType, VECTOR_LENGTH> cabana_particles; ///< Particle storage on device
  private:
    std::shared_ptr<Settings> settingsPtr;  ///< Pointer to Settings object
    ////////////////////////////////////////////////////////////////////////////
    // Cabana data structures (used for parallelization)
    // These are private and because data should be moved from device to host
    // before being accessed by the user. This can be implemented in public
    // methods, if necessarry
    Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>* simd_policy;                 ///< Policy used to access the particle data
    ListType neighbour_list; ///< Neighbour list
    Cabana::LinkedCellList<DeviceType> grid; ///< Grid used to accelerate the search for neighbors
    Kokkos::RangePolicy<ExecutionSpace> range_policy; ///< Policy used to loop over neighbors
    //using SerialHost = Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>;
    //Cabana::AoSoA<ParticleType, SerialHost, 8> particles_h("particles_h",n_particles); ///Temporary storage on host
    ////////////////////////////////////////////////////////////////////////////


    // used to track when a particle fails to find a solution in any EoS
    std::vector<int> particles_out_of_grid;

    std::vector<int> list;

    // eccentricities
    std::vector<double> timesteps;
    std::vector<double> e_2_X, e_2_P;

  public:
    // initialize system state, linklist, etc.
    void add_particle( Particle<D> p );
    void initialize();
    void allocate_cabana_particles();
    void initialize_linklist();
    void reset_neighbour_list();

    // check conserved quantities
    void conservation_entropy();
    void conservation_energy();
    void conservation_BSQ();

    //TODO: observables should probably be computed in a separate class
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

template<unsigned int D>
inline void SystemState<D>::add_particle( Particle<D> p )
{
  particles.push_back(p);
}
} // namespace ccake
#endif