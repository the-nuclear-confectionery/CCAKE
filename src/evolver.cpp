#include "evolver.h"

using namespace ccake;

//Template instantiations
template class Evolver<1>;
template class Evolver<2>;
template class Evolver<3>;

/// @brief Create object responsible for controlling hydro evolution.
/// @details This constructor sets up the evolver object with pointers to
/// the Settings and SystemState objects.
/// @tparam D The number of spatial dimensions.
/// @param settingsPtr_in Settings smart pointer.
/// @param systemPtr_in SystemState smart pointer.
template <unsigned int D>
Evolver<D>::Evolver(std::shared_ptr<Settings> settingsPtr_in,
                    std::shared_ptr<SystemState<D>> systemPtr_in)
                 : settingsPtr(settingsPtr_in),
                   systemPtr(systemPtr_in){};

/// @brief Allocate memory for the cache used for the Evolver.
/// @details The evolver needs to store a small amount of information of
/// the previous state of the system. This function allocates memory for
/// this cache.
/// @tparam D The number of spatial dimensions.
template <unsigned int D>
void Evolver<D>::allocate_cache()
{
    n_particles = systemPtr->n_particles;

    //Allocate on device memory for the system state before
    //evolution begins.
    evolver_cache = Cabana::AoSoA<EvolverCache, DeviceType, VECTOR_LENGTH>("cache", n_particles);
}

template <unsigned int D>
void Evolver<D>::allocate_k_values()
{
  n_particles = systemPtr->n_particles;
  k = Cabana::AoSoA<EvolverCache, DeviceType, VECTOR_LENGTH>("k", n_particles);
//   k2 = Cabana::AoSoA<EvolverCache, DeviceType, VECTOR_LENGTH>("k1", n_particles);
//   k3 = Cabana::AoSoA<EvolverCache, DeviceType, VECTOR_LENGTH>("k1", n_particles);
//   k4 = Cabana::AoSoA<EvolverCache, DeviceType, VECTOR_LENGTH>("k1", n_particles);
  auto kshv = Cabana::slice<evolver_cache_info::viscous_shear>(k);
  auto ku = Cabana::slice<evolver_cache_info::four_velocity>(k);
  auto kr = Cabana::slice<evolver_cache_info::position>(k);
  auto ks = Cabana::slice<evolver_cache_info::specific_entropy>(k);
  auto kBulk = Cabana::slice<evolver_cache_info::Bulk_pressure>(k);
  auto kE0 = Cabana::slice<evolver_cache_info::E0>(k);

  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());

  auto fill_k_cache = KOKKOS_LAMBDA (const int is, const int ia)
    {
      for (int i=0; i<D+1; i++)
      for (int j=0; j<D+1; j++)
        kshv.access(is, ia,i,j) = 0;
      for(int idir = 0; idir < D; ++idir)
      {
        ku.access(is, ia, idir) = 0;
        kr.access(is, ia, idir) = 0;
      }
      ks.access(is, ia) = 0;
      kBulk.access(is, ia) = 0;
      kE0.access(is, ia) = 0;
    };
    Cabana::simd_parallel_for(simd_policy, fill_k_cache, "update_k_cache");
    Kokkos::fence();
}

/// @brief Function to decide which algorith to use for the evolution.
/// @details This function is responsible for deciding which algorithm to use
/// for the evolution of the hydrodynamic system.
/// @note At the moment, only the Runge-Kutta 2nd order algorithm is
/// implemented. Additional algorithms could be added here.
/// @todo Instead of using a switch case, a map could be used to store the
/// algorithms and their names. This would allow for easier addition of new
/// algorithms.
/// @todo Instead of relying in the rk_order argument, we should use the
/// settingsPtr to decide which algorithm to use.
/// @tparam D The number of spatial dimensions.
/// @param dt
/// @param rk_order
/// @param time_derivatives_functional
template <unsigned int D>
void Evolver<D>::execute_timestep(double dt, int rk_order,
                                  std::function<void(void)> time_derivatives_functional )
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
    return;

}

/// @brief Fills the cache with the quantities at the current time step.
/// @details This function is responsible for filling the cache with the
/// quantities at the current time step, before we modify them while we
/// evolve the system.
/// @tparam D The number of spatial dimensions.
template <unsigned int D>
void Evolver<D>::set_current_timestep_quantities()
{

  CREATE_VIEW(device_, systemPtr->cabana_particles);
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());

  auto shv = Cabana::slice<evolver_cache_info::viscous_shear>(evolver_cache);
  auto u = Cabana::slice<evolver_cache_info::four_velocity>(evolver_cache);
  auto r = Cabana::slice<evolver_cache_info::position>(evolver_cache);
  auto s = Cabana::slice<evolver_cache_info::specific_entropy>(evolver_cache);
  auto Bulk = Cabana::slice<evolver_cache_info::Bulk_pressure>(evolver_cache);
  auto E0 = Cabana::slice<evolver_cache_info::E0>(evolver_cache);

  auto fill_cache = KOKKOS_LAMBDA (const int is, const int ia)
  {
    for (int i=0; i<D+1; i++)
    for (int j=0; j<D+1; j++)
      shv.access(is, ia,i,j) = device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, i, j);

    for(int idir = 0; idir < D; ++idir)
    {
      u.access(is, ia, idir) = device_hydro_vector.access(is, ia, ccake::hydro_info::u,idir);
      r.access(is, ia, idir) = device_position.access(is, ia, idir);
    }
    s.access(is, ia) = device_specific_density.access(is, ia, densities_info::s);
    Bulk.access(is, ia) = device_hydro_scalar.access(is, ia, hydro_info::Bulk);
    E0.access(is, ia) = device_contribution_to_total_Ez.access(is, ia);
  };
  Cabana::simd_parallel_for(simd_policy, fill_cache, "fill_cache");
}

///Calculates (k1 + 2*k2 + 2*k3 + k4) / 6
template <unsigned int D>
void Evolver<D>::update_k(int n)
{
  std::map<int, double> w = {{1, 1/6.0}, {2, 2/6.0}, {3, 2/6.0}, {4, 1/6.0}};
  CREATE_VIEW(device_, systemPtr->cabana_particles);
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());

  auto dshv_dt = Cabana::slice<evolver_cache_info::viscous_shear>(k);
  auto du_dt = Cabana::slice<evolver_cache_info::four_velocity>(k);
  auto dr_dt = Cabana::slice<evolver_cache_info::position>(k);
  auto ds_dt = Cabana::slice<evolver_cache_info::specific_entropy>(k);
  auto dBulk_dt = Cabana::slice<evolver_cache_info::Bulk_pressure>(k);
  auto dE0_dt = Cabana::slice<evolver_cache_info::E0>(k);

  auto update_k_cache = KOKKOS_LAMBDA (const int is, const int ia)
  {
    for (int i=0; i<D+1; i++)
    for (int j=0; j<D+1; j++)
      Kokkos::atomic_add(&dshv_dt.access(is, ia,i,j), w.at(n)*device_hydro_spacetime_matrix.access(is, ia, hydro_info::dshv_dt, i, j));

    for(int idir = 0; idir < D; ++idir)
    {
      Kokkos::atomic_add(&du_dt.access(is, ia, idir), w.at(n)*device_hydro_vector.access(is, ia, ccake::hydro_info::du_dt,idir));
      Kokkos::atomic_add(&dr_dt.access(is, ia, idir), w.at(n)*device_hydro_vector.access(is, ia, ccake::hydro_info::v, idir));
    }
    Kokkos::atomic_add(&ds_dt.access(is, ia), w.at(n)*device_d_dt_spec.access(is, ia, densities_info::s));;
    Kokkos::atomic_add(&dBulk_dt.access(is, ia), w.at(n)*device_hydro_scalar.access(is, ia, hydro_info::dBulk_dt));
    Kokkos::atomic_add(&dE0_dt.access(is, ia), w.at(n)*device_contribution_to_total_Ez.access(is, ia));
  };
  Cabana::simd_parallel_for(simd_policy, update_k_cache, "update_k_cache");
  Kokkos::fence();
}

///implements final updates for rk4 (y = y0 + (k1 + 2*k2 + 2*k3 + k4) * dt / 6)
template <unsigned int D>
void Evolver<D>::update_rk4(double dt){
      double w = 1/6.0;
  //Create views for the device
      CREATE_VIEW(device_, systemPtr->cabana_particles);
      auto slice_shv0 = Cabana::slice<evolver_cache_info::viscous_shear>(evolver_cache);
      auto slice_u0 = Cabana::slice<evolver_cache_info::four_velocity>(evolver_cache);
      auto slice_r0 = Cabana::slice<evolver_cache_info::position>(evolver_cache);
      auto slice_s0 = Cabana::slice<evolver_cache_info::specific_entropy>(evolver_cache);
      auto slice_Bulk0 = Cabana::slice<evolver_cache_info::Bulk_pressure>(evolver_cache);
      auto slice_E0 = Cabana::slice<evolver_cache_info::E0>(evolver_cache);
      //k<var> corresponds to w1k1 + 2*w2k2 + 2*w3k3 + w4k4 in RK4
      auto slice_kshv = Cabana::slice<evolver_cache_info::viscous_shear>(k);
      auto slice_ku = Cabana::slice<evolver_cache_info::four_velocity>(k);
      auto slice_kr = Cabana::slice<evolver_cache_info::position>(k);
      auto slice_ks = Cabana::slice<evolver_cache_info::specific_entropy>(k);
      auto slice_kBulk = Cabana::slice<evolver_cache_info::Bulk_pressure>(k);
      auto slice_kE0 = Cabana::slice<evolver_cache_info::E0>(k);


      //Add everything together
      auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());
      auto update_values = KOKKOS_LAMBDA(const int is, const int ia){
        double specific_s0 = slice_s0.access(is, ia);
        double ks = slice_ks.access(is, ia);
        double Bulk0 = slice_Bulk0.access(is, ia);
        double kBulk = slice_kBulk.access(is, ia);
        double particles_E0 = slice_E0.access(is, ia);
        double kE0 = slice_kE0.access(is, ia);

        milne::Matrix<double, D+1, D+1> shv0, kshv;
        milne::Vector<double, D> r0, kr, u0, ku;
        for(int idir=0; idir<D; ++idir){
          r0(idir) = slice_r0.access(is, ia, idir);
          kr(idir) = slice_kr.access(is, ia, idir);
          u0(idir) = slice_u0.access(is, ia, idir);
          ku(idir) = slice_ku.access(is, ia, idir);
        }
        for(int idir=1; idir<D+1; ++idir){
          for(int jdir=1; jdir<D+1; ++jdir){
            shv0(idir, jdir) = slice_shv0.access(is, ia, idir, jdir);
            kshv(idir, jdir) = slice_kshv.access(is, ia, idir, jdir);
          }
        }
        
        //Sum everything up
        double specific_s         = specific_s0         + dt*ks;
        double Bulk               = Bulk0               + dt*kBulk;
        double particles_Ez       = particles_E0        + dt*kE0;
        milne::Vector<double,D> r = r0                  + dt*kr;
        milne::Vector<double,D> u = u0                  + dt*ku;
        milne::Matrix<double,D,D> shv;
        for(int idir=0; idir<D; ++idir)
        for(int jdir=0; jdir<D; ++jdir)
          shv(idir,jdir)          = shv0(idir+1,jdir+1) + dt*kshv(idir, jdir);

        //Update in memory
        int freeze = device_freeze.access(is, ia); //Check if the cell is frozen
        if (specific_s < 0.0 && freeze > 3){ //If frozen, we do not want to crash because of negative entropy
          specific_s = 1.e-3; //Enforce positivity
        } else if (specific_s < 0.0){ //Else, something went terribly wrong
          formatted_output::detail("Negative entropy density");
          exit(EXIT_FAILURE);
        }

        device_specific_density.access(is, ia, densities_info::s) = specific_s;
        device_hydro_scalar.access(is, ia, hydro_info::Bulk) = Bulk;
        device_contribution_to_total_Ez.access(is, ia) = particles_Ez;
        for (int idir=0; idir<D; ++idir){
          device_position.access(is, ia, idir) = r(idir);
          device_hydro_vector.access(is, ia, hydro_info::u, idir) = u(idir);
        }
        for (int idir=1; idir<D+1; ++idir)
        for (int jdir=1; jdir<D+1; ++jdir)
            device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir) = shv(idir-1,jdir-1);

        //Enforce zero values for components greater than dimension D
        for (int idir=D; idir<3; ++idir){
          Kokkos::atomic_assign(&device_position.access(is, ia, idir) , 0.0);
          Kokkos::atomic_assign(&device_hydro_vector.access(is, ia, hydro_info::u, idir) , 0.0);
          for (int jdir=0; jdir<4; ++jdir){
            Kokkos::atomic_assign(&device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir+1, jdir), 0.0);
            Kokkos::atomic_assign(&device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, jdir, idir+1), 0.0);
          }
        }
      };
      Cabana::simd_parallel_for(simd_policy, update_values, "update_values");
      Kokkos::fence();

}

//==========================================================================
template <unsigned int D>
void Evolver<D>::step_rk(double dt, double t0, std::function<void(void)> time_derivatives_functional ){

  time_derivatives_functional();
  auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());

  //Create views for the device
  CREATE_VIEW(device_, systemPtr->cabana_particles);
  auto slice_shv0 = Cabana::slice<evolver_cache_info::viscous_shear>(evolver_cache);
  auto slice_u0 = Cabana::slice<evolver_cache_info::four_velocity>(evolver_cache);
  auto slice_r0 = Cabana::slice<evolver_cache_info::position>(evolver_cache);
  auto slice_s0 = Cabana::slice<evolver_cache_info::specific_entropy>(evolver_cache);
  auto slice_Bulk0 = Cabana::slice<evolver_cache_info::Bulk_pressure>(evolver_cache);
  auto slice_E0 = Cabana::slice<evolver_cache_info::E0>(evolver_cache);

  auto update_rk2_step = KOKKOS_LAMBDA(const int is, const int ia){
    //Cache previous step values
    double specific_s0 = slice_s0.access(is, ia);
    double Bulk0 = slice_Bulk0.access(is, ia);
    double particles_E0 = slice_E0.access(is, ia);
    milne::Matrix<double, D+1, D+1> shv0;
    milne::Vector<double, D> r0, u0;
    for(int idir=0; idir<D; ++idir){
      r0(idir) = slice_r0.access(is, ia, idir);
      u0(idir) = slice_u0.access(is, ia, idir);
    }
    for(int idir=1; idir<D+1; ++idir)
    for(int jdir=1; jdir<D+1; ++jdir)
      shv0(idir, jdir) = slice_shv0.access(is, ia, idir, jdir);

    //Cache derivatives
    double dBulk_dt = device_hydro_scalar.access(is,ia,hydro_info::dBulk_dt);
    double d_dt_specific_s = device_d_dt_spec.access(is, ia, densities_info::s);
    double dEz_dt = device_contribution_to_total_dEz.access(is,ia);
    milne::Vector<double, D> du_dt, v;
    for(int idir=0; idir<D; ++idir){
      du_dt(idir) = device_hydro_vector.access(is, ia, hydro_info::du_dt, idir);
      v(idir) = device_hydro_vector.access(is, ia, hydro_info::v, idir);
    }
    milne::Matrix<double, D, D> dshv_dt;
    for(int idir=0; idir<D; ++idir)
    for(int jdir=0; jdir<D; ++jdir)
      dshv_dt(idir, jdir) = device_hydro_space_matrix.access(is, ia, hydro_info::dshv_dt, idir, jdir);

    //Compute updated quantities
    double specific_s         = specific_s0         + dt*d_dt_specific_s;
    double Bulk               = Bulk0               + dt*dBulk_dt;
    double particles_Ez       = particles_E0        + dt*dEz_dt;
    milne::Vector<double,D> r = r0                  + dt*v;
    milne::Vector<double,D> u = u0                  + dt*du_dt;
    milne::Matrix<double,D,D> shv;
    for(int idir=0; idir<D; ++idir)
    for(int jdir=0; jdir<D; ++jdir)
      shv(idir,jdir)          = shv0(idir+1,jdir+1) + dt*dshv_dt(idir, jdir);

    //Update in memory
    int freeze = device_freeze.access(is, ia); //Check if the cell is frozen
    if (specific_s < 0.0 && freeze > 3){ //If frozen, we do not want to crash because of negative entropy
      specific_s = 1.e-3; //Enforce positivity
    } else if (specific_s < 0.0){ //Else, something went terribly wrong
      formatted_output::detail("Negative entropy density");
      exit(EXIT_FAILURE);
    }

    device_specific_density.access(is, ia, densities_info::s) = specific_s;
    device_hydro_scalar.access(is, ia, hydro_info::Bulk) = Bulk;
    device_contribution_to_total_Ez.access(is, ia) = particles_Ez;
    for (int idir=0; idir<D; ++idir){
      device_position.access(is, ia, idir) = r(idir);
      device_hydro_vector.access(is, ia, hydro_info::u, idir) = u(idir);
    }
    for (int idir=1; idir<D+1; ++idir)
    for (int jdir=1; jdir<D+1; ++jdir)
        device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir) = shv(idir-1,jdir-1);

    //Enforce zero values for components greater than dimension D
    for (int idir=D; idir<3; ++idir){
      Kokkos::atomic_assign(&device_position.access(is, ia, idir) , 0.0);
      Kokkos::atomic_assign(&device_hydro_vector.access(is, ia, hydro_info::u, idir) , 0.0);
      for (int jdir=0; jdir<4; ++jdir){
        Kokkos::atomic_assign(&device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir+1, jdir), 0.0);
        Kokkos::atomic_assign(&device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, jdir, idir+1), 0.0);
      }
    }
  };
  Cabana::simd_parallel_for(simd_policy, update_rk2_step, "update_rk2_step");
  Kokkos::fence();
  systemPtr->t  = t0 + dt;
  double t = systemPtr->t;
  Cabana::simd_parallel_for(simd_policy, KOKKOS_LAMBDA(const int is, const int ia)
  {
    device_hydro_scalar.access(is, ia,ccake::hydro_info::t) = t;
  }, "update_hydro_time_step");
  Kokkos::fence();
}

/// @brief Evolve the system using a second order Runge-Kutta algorithm.
/// @tparam D The number of spatial dimensions.
/// @param dt The time step to be used in the evolution.
/// @param time_derivatives_functional A functional that computes the time
/// derivatives of the system.
template <unsigned int D>
void Evolver<D>::advance_timestep_rk2( double dt,
                                        std::function<void(void)> time_derivatives_functional )
{
      // initialize quantities at current time step
      set_current_timestep_quantities();
      double t0      = systemPtr->t;
      ////////////////////////////////////////////
      //    first step
      ////////////////////////////////////////////
      formatted_output::report("RK(n=2) evolution, step 1");
      step_rk(.5*dt, t0, time_derivatives_functional);
      // E1   = dt*systemPtr->dEz;

      ////////////////////////////////////////////
      //    second step
      ////////////////////////////////////////////
      formatted_output::report("RK(n=2) evolution, step 2");
      step_rk(dt, t0, time_derivatives_functional);
      // E2   = dt*systemPtr->dEz;
      // constexpr double w1 = 1.0/6.0, w2 = 1.0/3.0;
      // systemPtr->Ez = E0 + w1*E1 + w2*E2;

      return;
}

template <unsigned int D>
void Evolver<D>::advance_timestep_rk4( double dt,
                                        std::function<void(void)> time_derivatives_functional )
{
      // initialize quantities at current time step
      set_current_timestep_quantities();
      allocate_k_values();
      double t0      = systemPtr->t;
      ////////////////////////////////////////////
      //    first step
      ////////////////////////////////////////////
      formatted_output::report("RK(n=4) evolution, step 1");
      step_rk(.5*dt, t0, time_derivatives_functional);
      

      ////////////////////////////////////////////
      //    second step
      ////////////////////////////////////////////
      formatted_output::report("RK(n=4) evolution, step 2");
      step_rk(0.5*dt, t0, time_derivatives_functional);
      update_k(1);

      ////////////////////////////////////////////
      //    third step
      ////////////////////////////////////////////
      formatted_output::report("RK(n=4) evolution, step 3");
      step_rk(0.5*dt, t0, time_derivatives_functional);
      update_k(2);

      ////////////////////////////////////////////
      //    fourth step
      ////////////////////////////////////////////
      formatted_output::report("RK(n=4) evolution, step 4");
      step_rk(dt, t0, time_derivatives_functional);
      update_k(3);

      time_derivatives_functional();
      update_k(3);

      update_rk4(dt);

      return;
}

