#ifndef EVOLVER_H
#define EVOLVER_H

#include <algorithm>
#include <string>
#include <vector>

#include "eos.h"
#include "formatted_output.h"
#include "particle.h"
#include "system_state.h"

using std::string;

class Evolver
{
  private:

    static constexpr bool REGULATE_LOW_T = false;

    
    Settings    * settingsPtr = nullptr;
    SystemState * systemPtr   = nullptr;

    int n_particles = -1;

    // creating vectors of quantities for RK evolution
    vector<double> specific_s0;
    vector<double> Bulk0;
    vector<double> particles_E0;

    vector< Vector<double,2> > u0;
    vector< Vector<double,2> > r0;

    vector< Matrix <double,2,2> > shv0;

  public:

    void set_SettingsPtr(Settings * settingsPtr_in) { settingsPtr = settingsPtr_in; }
    void set_SystemStatePtr(SystemState * systemPtr_in) { systemPtr = systemPtr_in; }

    //==========================================================================
    void advance_timestep( double dt, int rk_order,
                            std::function<void(void)>
                            time_derivatives_functional )
    {
      switch ( rk_order )
      {
        case 2:
          advance_timestep_rk2( dt, time_derivatives_functional );
          break;
        case 4:
          advance_timestep_rk4( dt, time_derivatives_functional );
          break;
        default:
          std::cerr << "Invalid Runge-Kutta order!" << std::endl;
          exit(8);
          break;
      }
    }




    //==============================================================================
    // this routine is used to initialize quantities prior to RK evolution
    void set_current_timestep_quantities()
    {
      // set number of particles
      n_particles = systemPtr->particles.size();

      // store quantities at current timestep
      specific_s0.resize(n_particles);
      Bulk0.resize(n_particles);
      particles_E0.resize(n_particles);
      u0.resize(n_particles);
      r0.resize(n_particles);
      shv0.resize(n_particles);

      for (int i = 0; i < n_particles; ++i)
      {
        auto & p = systemPtr->particles[i];

        u0[i]          = p.hydro.u;
        r0[i]          = p.r;
        specific_s0[i] = p.specific.s;
        Bulk0[i]       = p.hydro.Bulk;
        mini( shv0[i], p.hydro.shv );

        particles_E0[i] = p.contribution_to_total_Ez;
      }
    }




    //==========================================================================
    void advance_timestep_rk2( double dt, std::function<void(void)>
                                          time_derivatives_functional )
    {
      systemPtr->rk2 = 1;
      double t0      = systemPtr->t;
      double E0      = systemPtr->Ez;

      // initialize quantities at current time step
      set_current_timestep_quantities();

      ////////////////////////////////////////////
      //    first step
      ////////////////////////////////////////////
      formatted_output::report("RK(n=2) evolution, step 1");

      // compute derivatives
      time_derivatives_functional();

      // update quantities
      {
        for (int i = 0; i < n_particles; i++)
        {
          auto & p    = systemPtr->particles[i];
          auto & ph   = p.hydro;

          p.r = r0[i] + 0.5*dt*ph.v;
          ph.u            = u0[i]          + 0.5*dt*ph.du_dt;
          p.specific.s    = specific_s0[i] + 0.5*dt*p.d_dt_spec.s;
          ph.Bulk         = Bulk0[i]       + 0.5*dt*ph.dBulk_dt;
          tmini( ph.shv,    shv0[i]        + 0.5*dt*ph.dshv_dt );


          // regulate updated results if necessary
          if ( REGULATE_LOW_T && p.specific.s < 0.0
                /*&& p.T() < 50.0/constants::hbarc_MeVfm*/ )
          {
            std::cerr << "WARNING: regulating particle " << p.ID << " at e = "
                      << p.e()*constants::hbarc_MeVfm << " MeV/fm^3 "
                      << "(T = " << p.T()*constants::hbarc_MeVfm << " MeV): "
                      << p.specific.s << " --> " << 0.5*(specific_s0[i]+1.0) << endl;
            p.specific.s    = 0.5*(specific_s0[i]+1.0);
          }


          p.contribution_to_total_Ez = particles_E0[i]
                                        + 0.5*dt*p.contribution_to_total_dEz;
        }
      }

      systemPtr->Ez = E0 + 0.5*dt*systemPtr->dEz;
      systemPtr->t  = t0 + 0.5*dt;

      ////////////////////////////////////////////
      //    second step
      ////////////////////////////////////////////
      formatted_output::report("RK(n=2) evolution, step 2");

      // compute derivatives
      time_derivatives_functional();

      // update quantities
      {
        for (int i = 0; i < n_particles; i++)
        {
          auto & p    = systemPtr->particles[i];
          auto & ph   = p.hydro;

          p.r = r0[i] + dt*ph.v;
            ph.u            = u0[i]          + dt*ph.du_dt;
            p.specific.s    = specific_s0[i] + dt*p.d_dt_spec.s;
            ph.Bulk         = Bulk0[i]       + dt*ph.dBulk_dt;
            tmini( ph.shv,    shv0[i]        + dt*ph.dshv_dt );

          // regulate updated results if necessary
          if ( REGULATE_LOW_T && p.specific.s < 0.0
                /*&& p.T() < 50.0/constants::hbarc_MeVfm*/ )
          {
            std::cerr << "WARNING: regulating particle " << p.ID << " at e = "
                      << p.e()*constants::hbarc_MeVfm << " MeV/fm^3 "
                      << "(T = " << p.T()*constants::hbarc_MeVfm << " MeV): "
                      << p.specific.s << " --> " << 0.5*(specific_s0[i]+1.0) << endl;
            p.specific.s    = 0.5*(specific_s0[i]+1.0);
          }


          p.contribution_to_total_Ez = particles_E0[i]
                                        + dt*p.contribution_to_total_dEz;

    //cout << "CHECK energies: " << i << "   " << t0+dt << "   "
    //      << p.r << "   " << p.e() << "   "
    //      << particles_E0[i] << "   "
    //      << p.contribution_to_total_E << endl;
        }
      }

      systemPtr->Ez = E0 + dt*systemPtr->dEz;
      systemPtr->t  = t0 + dt;

      return;
    }



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
        ps.ets1      = dt*p.d_dt_spec.s;
        ps.b1        = dt*ph.dBulk_dt;
        ps.shv1      = dt*ph.dshv_dt;

        // implement increments with appropriate coefficients
        ph.u         = u0[i]          + 0.5*ps.k1;
        p.r          = r0[i]          + 0.5*ps.r1;
        p.specific.s = specific_s0[i] + 0.5*ps.ets1;
        ph.Bulk      = Bulk0[i]       + 0.5*ps.b1;
        tmini(ph.shv,  shv0[i]        + 0.5*ps.shv1);

        // regulate updated results if necessary
        if ( REGULATE_LOW_T && p.specific.s < 0.0
              && p.T() < 50.0/constants::hbarc_MeVfm )
          p.specific.s    = specific_s0[i];
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
        ps.ets2      = dt*p.d_dt_spec.s;
        ps.b2        = dt*ph.dBulk_dt;
        ps.shv2      = dt*ph.dshv_dt;

        ph.u         = u0[i]          + 0.5*ps.k2;
        p.r          = r0[i]          + 0.5*ps.r2;
        p.specific.s = specific_s0[i] + 0.5*ps.ets2;
        ph.Bulk      = Bulk0[i]       + 0.5*ps.b2;
        tmini(ph.shv,  shv0[i]        + 0.5*ps.shv2);

        // regulate updated results if necessary
        if ( REGULATE_LOW_T && p.specific.s < 0.0
              && p.T() < 50.0/constants::hbarc_MeVfm )
          p.specific.s    = specific_s0[i];
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
        ps.ets3      = dt*p.d_dt_spec.s;
        ps.b3        = dt*ph.dBulk_dt;
        ps.shv3      = dt*ph.dshv_dt;


        ph.u         = u0[i]          + ps.k3;
        p.r          = r0[i]          + ps.r3;
        p.specific.s = specific_s0[i] + ps.ets3;
        ph.Bulk      = Bulk0[i]       + ps.b3;
        tmini(ph.shv,  shv0[i]        + ps.shv3);

        // regulate updated results if necessary
        if ( REGULATE_LOW_T && p.specific.s < 0.0
              && p.T() < 50.0/constants::hbarc_MeVfm )
          p.specific.s    = specific_s0[i];
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
        ps.ets4      = dt*p.d_dt_spec.s;
        ps.b4        = dt*ph.dBulk_dt;
        ps.shv4      = dt*ph.dshv_dt;

        // sum the weighted steps into yf and return the final y values
        ph.u         = u0[i]          + w1*ps.k1   + w2*ps.k2   + w2*ps.k3   + w1*ps.k4;
        p.r          = r0[i]          + w1*ps.r1   + w2*ps.r2   + w2*ps.r3   + w1*ps.r4;
        p.specific.s = specific_s0[i] + w1*ps.ets1 + w2*ps.ets2 + w2*ps.ets3 + w1*ps.ets4;
        ph.Bulk      = Bulk0[i]       + w1*ps.b1   + w2*ps.b2   + w2*ps.b3   + w1*ps.b4;
        tmini(ph.shv,  shv0[i]        + w1*ps.shv1 + w2*ps.shv2 + w2*ps.shv3 + w1*ps.shv4);

        // regulate updated results if necessary
        if ( REGULATE_LOW_T && p.specific.s < 0.0
              && p.T() < 50.0/constants::hbarc_MeVfm )
          p.specific.s    = specific_s0[i];
      }

      E4            = dt*systemPtr->dEz;
      systemPtr->Ez = E0 + w1*E1 + w2*E2 + w2*E3 + w1*E4;


    }





};

typedef std::shared_ptr<Evolver> pEvolver;  // smart pointer to evolver object

#endif