#ifndef EVOLVER_H
#define EVOLVER_H

#include <algorithm>
#include <string>
#include <vector>

#include <Cabana_Core.hpp>

#include "eos.h"
#include "formatted_output.h"
#include "particle.h"
#include "system_state.h"
#include "utilities.h"
#include "milne.hpp"
#include "cartesian.hpp"

namespace ccake{


using EvolverCache = Cabana::MemberTypes<double[4][4], // stress-energy tensor
                                         double[3][4], // diffusion 
                                         double[3],    // four-velocity
                                         double[3],    // position
                                         double,       // extensive_entropy
                                         double,       // extensive_baryon
                                         double,       // extensive_strangeness
                                         double,       // extensive_electric
                                         double,       // Bulk pressure
                                         double        // E0 = Ez
>;

namespace evolver_cache_info
{
enum cache_data
{
  extensive_shear,
  extensive_diffusion,
  four_velocity,
  position,
  extensive_entropy,
  extensive_baryon,
  extensive_strangeness,
  extensive_electric,
  extensive_bulk,
  E0
};
}

/// @class Evolver
/// @brief Class responsible for evolving the hydrodynamic system in time
/// @details This class is responsible for evolving the hydrodynamic system
/// in time. Variables to be evolved are extensive entropy, space components
/// of four-velocity, position, bulk pressure, and the transverse components
/// of the viscous shear tensor.
/// @tparam D
template <unsigned int D>
class Evolver
{
  private:

    static constexpr bool REGULATE_NEGATIVE_S = false;
    static constexpr bool REGULATE_LARGE_S    = false;


    std::shared_ptr<Settings> settingsPtr;
    std::shared_ptr<SystemState<D>> systemPtr;

    int n_particles = -1;

    Cabana::AoSoA<EvolverCache, DeviceType, VECTOR_LENGTH> evolver_cache;
    Cabana::AoSoA<EvolverCache, DeviceType, VECTOR_LENGTH> k;


    // creating vectors of quantities for RK evolution
    vector<double> extensive_s0;
    vector<double> extensive_rhoB0;
    vector<double> extensive_rhoS0;
    vector<double> extensive_rhoQ0;

    vector<double> rhoB0;
    vector<double> rhoS0;
    vector<double> rhoQ0;
    
    vector<double> Bulk0;
    vector<double> particles_E0;

    vector< Vector<double,D> > u0;
    vector< Vector<double,D> > r0;
    vector< Matrix<double,3,4> > diffusion0;
    vector< Matrix <double,4,4> > shv0;

    void advance_timestep_rk2( double dt,
                               std::function<void(void)> time_derivatives_functional );
    void advance_timestep_rk4( double dt,
                               std::function<void(void)> time_derivatives_functional );


  public:

    Evolver() = delete;
    Evolver(std::shared_ptr<Settings> settingsPtr_in,
            std::shared_ptr<SystemState<D>> systemPtr_in);


    //==========================================================================
    void execute_timestep(double dt, int rk_order,
                          std::function<void(void)> time_derivatives_functional );
    void set_current_timestep_quantities();
    void update_k(int n);
    void allocate_cache();
    void allocate_k_values();
    void step_rk(double dt, double t0, std::function<void(void)> time_derivatives_functional );
    void update_rk4(double dt);

};
}
#endif

/*

    //==========================================================================
    void advance_timestep_rk4( double dt, std::function<void(void)>
                                          time_derivatives_functional )
    {
      // define a local struct just to help with the RK evolution here
      struct particle_state_RK4
      {
        double ets1 = 0.0, ets2 = 0.0, ets3 = 0.0, ets4 = 0.0;
        double b1 = 0.0, b2 = 0.0, b3 = 0.0, b4 = 0.0;
        Vector<double,2> k1, k2, k3, k4;
        Vector<double,2> r1, r2, r3, r4;
        Matrix<double,2,2> shv1, shv2, shv3, shv4;
      };

      // make vector of struct objects
      const int number_of_particles = systemPtr->particles.size();
      vector<particle_state_RK4> particle_states( number_of_particles );

      // set up RK4 evolution
      systemPtr->rk2 = 1;
      double t0      = systemPtr->t;
      double E0      = systemPtr->Ez;

      double E1 = 0.0, E2 = 0.0, E3 = 0.0, E4 = 0.0;

      // initialize quantities at current time step
      set_current_timestep_quantities();

      ////////////////////////////////////////////
      //    first step
      ////////////////////////////////////////////
      formatted_output::report("RK(n=4) evolution, step 1");

      // compute derivatives
      time_derivatives_functional();

      for (int i = 0; i < number_of_particles; i++)
      {
        auto & p    = systemPtr->particles[i];
        auto & ph   = p.hydro;
        auto & ps   = particle_states[i];

        // store increments
        ps.k1        = dt*ph.du_dt;
        ps.r1        = dt*ph.v;
        ps.ets1      = dt*p.d_dt_extensive.s;
        ps.b1        = dt*ph.dBulk_dt;
        ps.shv1      = dt*ph.dshv_dt;

        // implement increments with appropriate coefficients
        ph.u         = u0[i]          + 0.5*ps.k1;
        p.r          = r0[i]          + 0.5*ps.r1;
        p.extensive.s = extensive_s0[i] + 0.5*ps.ets1;
        ph.Bulk      = Bulk0[i]       + 0.5*ps.b1;
        tmini(ph.shv,  shv0[i]        + 0.5*ps.shv1);

        // regulate updated results if necessary
        if ( REGULATE_NEGATIVE_S && p.extensive.s < 0.0
              && p.T() < 50.0/constants::hbarc_MeVfm )
          p.extensive.s    = extensive_s0[i];
      }

      E1           = dt*systemPtr->dEz;


      ////////////////////////////////////////////
      //    second step
      ////////////////////////////////////////////
      formatted_output::report("RK(n=4) evolution, step 2");

      systemPtr->t = t0 + 0.5*dt;
      time_derivatives_functional();

      for (int i = 0; i < number_of_particles; i++)
      {
        auto & p    = systemPtr->particles[i];
        auto & ph   = p.hydro;
        auto & ps   = particle_states[i];

        ps.k2        = dt*ph.du_dt;
        ps.r2        = dt*ph.v;
        ps.ets2      = dt*p.d_dt_extensive.s;
        ps.b2        = dt*ph.dBulk_dt;
        ps.shv2      = dt*ph.dshv_dt;

        ph.u         = u0[i]          + 0.5*ps.k2;
        p.r          = r0[i]          + 0.5*ps.r2;
        p.extensive.s = extensive_s0[i] + 0.5*ps.ets2;
        ph.Bulk      = Bulk0[i]       + 0.5*ps.b2;
        tmini(ph.shv,  shv0[i]        + 0.5*ps.shv2);

        // regulate updated results if necessary
        if ( REGULATE_NEGATIVE_S && p.extensive.s < 0.0
              && p.T() < 50.0/constants::hbarc_MeVfm )
          p.extensive.s    = extensive_s0[i];
      }

      E2           = dt*systemPtr->dEz;


      ////////////////////////////////////////////
      //    third step
      ////////////////////////////////////////////
      formatted_output::report("RK(n=4) evolution, step 3");

      time_derivatives_functional();

      for (int i = 0; i < number_of_particles; i++)
      {
        auto & p    = systemPtr->particles[i];
        auto & ph   = p.hydro;
        auto & ps   = particle_states[i];

        ps.k3        = dt*ph.du_dt;
        ps.r3        = dt*ph.v;
        ps.ets3      = dt*p.d_dt_extensive.s;
        ps.b3        = dt*ph.dBulk_dt;
        ps.shv3      = dt*ph.dshv_dt;


        ph.u         = u0[i]          + ps.k3;
        p.r          = r0[i]          + ps.r3;
        p.extensive.s = extensive_s0[i] + ps.ets3;
        ph.Bulk      = Bulk0[i]       + ps.b3;
        tmini(ph.shv,  shv0[i]        + ps.shv3);

        // regulate updated results if necessary
        if ( REGULATE_NEGATIVE_S && p.extensive.s < 0.0
              && p.T() < 50.0/constants::hbarc_MeVfm )
          p.extensive.s    = extensive_s0[i];
      }

      E3           = dt*systemPtr->dEz;

      ////////////////////////////////////////////
      //    fourth step
      ////////////////////////////////////////////
      formatted_output::report("RK(n=4) evolution, step 4");

      systemPtr->t = t0 + dt;
      time_derivatives_functional();

      constexpr double w1 = 1.0/6.0, w2 = 1.0/3.0;
      for (int i = 0; i < number_of_particles; i++)
      {
        auto & p    = systemPtr->particles[i];
        auto & ph   = p.hydro;
        auto & ps   = particle_states[i];

        ps.k4        = dt*ph.du_dt;
        ps.r4        = dt*ph.v;
        ps.ets4      = dt*p.d_dt_extensive.s;
        ps.b4        = dt*ph.dBulk_dt;
        ps.shv4      = dt*ph.dshv_dt;

        // sum the weighted steps into yf and return the final y values
        ph.u         = u0[i]          + w1*ps.k1   + w2*ps.k2   + w2*ps.k3   + w1*ps.k4;
        p.r          = r0[i]          + w1*ps.r1   + w2*ps.r2   + w2*ps.r3   + w1*ps.r4;
        p.extensive.s = extensive_s0[i] + w1*ps.ets1 + w2*ps.ets2 + w2*ps.ets3 + w1*ps.ets4;
        ph.Bulk      = Bulk0[i]       + w1*ps.b1   + w2*ps.b2   + w2*ps.b3   + w1*ps.b4;
        tmini(ph.shv,  shv0[i]        + w1*ps.shv1 + w2*ps.shv2 + w2*ps.shv3 + w1*ps.shv4);

        // regulate updated results if necessary
        if ( REGULATE_NEGATIVE_S && p.extensive.s < 0.0
              && p.T() < 50.0/constants::hbarc_MeVfm )
          p.extensive.s    = extensive_s0[i];
      }

      E4            = dt*systemPtr->dEz;
      systemPtr->Ez = E0 + w1*E1 + w2*E2 + w2*E3 + w1*E4;


    }

*/