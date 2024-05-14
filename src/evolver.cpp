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
          cout << "RK4 Not implemented" << std::endl;
          Kokkos::finalize();
          exit(-42);
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
  //auto E0 = Cabana::slice<evolver_cache_info::E0>(evolver_cache);
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
    //E0.access(is, ia) = device_contribution_to_total_E.access(is, ia);
  };
  Cabana::simd_parallel_for(simd_policy, fill_cache, "fill_cache");
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

  auto update_rk2_step = KOKKOS_LAMBDA(const int is, const int ia){
    //Cache previous step values
    double specific_s0 = slice_s0.access(is, ia);
    double Bulk0 = slice_Bulk0.access(is, ia);
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
    milne::Vector<double,D> r = r0                  + dt*v;
    milne::Vector<double,D> u = u0                  + dt*du_dt;
    milne::Matrix<double,D,D> shv;
    for(int idir=0; idir<D; ++idir)
    for(int jdir=0; jdir<D; ++jdir)
      shv(idir,jdir)          = shv0(idir+1,jdir+1) + dt*dshv_dt(idir, jdir);

    device_specific_density.access(is, ia, densities_info::s) = specific_s;
    device_hydro_scalar.access(is, ia, hydro_info::Bulk) = Bulk;
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

  if (settingsPtr->regulate_dissipative_terms ){
    auto regulate = KOKKOS_LAMBDA(const int is, const int ia){
      //Enforce minimum value for specific entropy
      int freeze = device_freeze.access(is, ia); //Check if the cell is frozen
      double specific_s = device_specific_density.access(is, ia, densities_info::s);
      if (specific_s < 0.0 && freeze > 0){ //If frozen, we do not want to crash because of negative entropy
        device_specific_density.access(is, ia, densities_info::s) = 1.e-3; //Enforce positivity
      } else if (specific_s < 0.0){
        exit(EXIT_FAILURE);
      }
    };
    Cabana::simd_parallel_for(simd_policy, regulate, "regulate");
    Kokkos::fence();
  }
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
      //double E0      = systemPtr->Ez;
      //Kokkos::View<int, MemorySpace> E0("E0");
      //Kokkos::View<int, MemorySpace> E("E");
      //Kokkos::deep_copy(E0, systemPtr->Ez);
      //Kokkos::deep_copy(E, systemPtr->Ez);
      // initialize quantities at current time step
      set_current_timestep_quantities();
      double t0      = systemPtr->t;
      ////////////////////////////////////////////
      //    first step
      ////////////////////////////////////////////
      step_rk(.5*dt, t0, time_derivatives_functional);

      ////////////////////////////////////////////
      //    second step
      ////////////////////////////////////////////
      step_rk(dt, t0, time_derivatives_functional);
      return;
}