#ifndef SYSTEM_STATE_H
#define SYSTEM_STATE_H

#include <Cabana_Core.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <array>

#include "constants.h"
#include "eos.h"
#include "formatted_output.h"
#include "kernel.h"
#include "particle.h"
#include "settings.h"
#include "stopwatch.h"
#include "vector.h"
#include "utilities.h"
#include "densities.h"
#include "milne.hpp"

using std::string;
using std::vector;


namespace ccake{
/// @class SystemState
/// @brief A class to hold the state of the system.
/// @details This class holds the state of the system, including the particles,
/// the equation of state, the settings, the time, the timestep, the total
/// entropy, energy, baryon charge, strangeness, electric charge, etc. It also
/// holds the Cabana data structures used for parallelization.
/// @todo: Computation of eccentricities and total energy needs to be implemented
/// in parallel (using Kokkos and Cabana)
template <unsigned int D>
class SystemState
{

  public:

    SystemState() = delete; ///< Default constructor is deleted. Settings must be passed in.
    SystemState( shared_ptr<Settings> settingsPtr_in): settingsPtr(settingsPtr_in) {};
     ~SystemState(){}

    void set_EquationOfStatePtr( EquationOfState * eosPtr_in );
    void set_SettingsPtr( Settings * settingsPtr_in );

    double t          = 0.0;

    double S          = 0.0; ///<Total entropy
    double S0         = 0.0; ///<Initial total entropy
    double Btotal     = 0.0; ///<Total baryon charge
    double Btotal0    = 0.0; ///<Initial total baryon charge
    double Stotal     = 0.0; ///<Total strangeness charge
    double Stotal0    = 0.0; ///<Initial total strangeness charge
    double Qtotal     = 0.0; ///<Total electric charge
    double Qtotal0    = 0.0; ///<Initial total electric charge

    double E          = 0.0;
    double Ez         = 0.0;
    double E0         = 0.0;
    double Eloss      = 0.0;
    double dEz        = 0.0;
    double Etot       = 0.0;
    double Ejet       = 0.0;
    double Ejet0      = 0.0;

    // double e_2_P      = 0.0;
    // double e_2_X      = 0.0;
    // int count_X       = 0;
    // int count_P       = 0;

    double efcheck    = 0.0;
    double sfcheck    = 0.0;

   /* INTS */
    int number_part_fo = 0;
    int n_particles = 0;

    int number_of_elapsed_timesteps = 0;


    std::vector<Particle<D>> particles;     ///< Vector of particles
    ////////////////////////////////////////////////////////////////////////////
    // Cabana data structures (used for parallelization)
    // These are private and because data should be moved from device to host
    // before being accessed by the user. This can be implemented in public
    // methods, if necessary
    //Cabana::LinkedCellList<DeviceType> grid; ///< Grid used to accelerate the search for neighbors //< I think this may not be necessary
    //using SerialHost = Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>;
    //Cabana::AoSoA<ParticleType, SerialHost, 8> particles_h("particles_h",n_particles); ///Temporary storage on host
    ////////////////////////////////////////////////////////////////////////////
    Cabana::AoSoA<CabanaParticle, DeviceType, VECTOR_LENGTH> cabana_particles; ///< Particle storage on device
    ListType neighbour_list; ///< Neighbour list
    std::vector<double> timesteps;
    std::vector<std::vector<double>> e_2_P_history_by_slice, e_2_X_history_by_slice;
    std::vector<std::vector<int>> count_P_history_by_slice, count_X_history_by_slice;
    std::vector<double> eta_slices;
  private:
    std::shared_ptr<Settings> settingsPtr;  ///< Pointer to Settings object

    // eccentricities
    // std::vector<double> timesteps;
    // // std::vector<double> e_2_X_history, e_2_P_history;
    // std::vector<std::vector<double>> e_2_P_history_by_slice, e_2_X_history_by_slice;
    // std::vector<std::vector<int>> count_P_history_by_slice, count_X_history_by_slice;
    // std::vector<double> eta_slices;
    // std::vector<int> count_X_history, count_P_history;


  public:
    // initialize system state, linklist, etc.
    void add_particle( Particle<D> p );
    void initialize();
    void allocate_cabana_particles();
    void copy_device_to_host();
    void copy_host_to_device();
    void initialize_linklist();
    void reset_neighbour_list();

    // check conserved quantities
    void conservation_entropy(bool first_iteration=false);
    void conservation_energy(bool first_iteration=false,double t = 1.0);
    void conservation_BSQ(bool first_iteration=false);

    //TODO: observables should probably be computed in a separate class
    void compute_eccentricities();
    // void compute_e_2_P(double slice);
    // void compute_e_2_X(double slice);
    void compute_e_2_P(double slice, int i);
    void compute_e_2_X(double slice, int i);

    int n(){ return n_particles; }
    double get_particle_T(int id) {return particles[id].T();}
    double get_particle_Freeze(int id) {return particles[id].Freeze;}
    void print_neighbors(int idx){
    //Loop of neighbors, showing them
      CREATE_VIEW(device_,cabana_particles)
        std::cout << "System has " << cabana_particles.size() << " particles" << std::endl;
        int num_n = Cabana::NeighborList<ListType>::numNeighbor( neighbour_list, idx );
	//if (num_n << 5) abort();
        std::cout << "Particle " << idx << " # neighbor = " << num_n << std::endl;
        for ( int j = 0; j < num_n; ++j ){
          int neighIdx = Cabana::NeighborList<ListType>::getNeighbor( neighbour_list, idx, j );
          std::cout << "    neighbor " << j << " = "
                    << neighIdx  << ", extensive.s: " << device_extensive(neighIdx, ccake::densities_info::s)
                    << ", sph_mass.s: " << device_sph_mass(neighIdx, ccake::densities_info::s) << std::endl;
                }
    };

    std::vector<std::array<double, 4>> get_particle_data(int idx){
      CREATE_VIEW(device_, cabana_particles);
      int num_neighbors = Cabana::NeighborList<ListType>::numNeighbor(neighbour_list, idx);
      std::vector<std::array<double, 4>> result;
      result.reserve(1);
      result.push_back({
          static_cast<double>(num_neighbors),         // number of neighbors
          device_position(idx, 0),                    // x
          device_position(idx, 1),                    // y
          device_position(idx, 2)                     // eta
      });

      return result;
    };

    // };



};

template<unsigned int D>
inline void SystemState<D>::add_particle( Particle<D> p )
{
  particles.push_back(p);
}
} // namespace ccake
#endif
