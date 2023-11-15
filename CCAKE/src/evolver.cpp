#include "evolver.h"

using namespace ccake;

//Template instantiations
template class Evolver<1>;
template class Evolver<2>;
template class Evolver<3>;

template <unsigned int D>
Evolver<D>::Evolver(std::shared_ptr<Settings> settingsPtr_in,
                    std::shared_ptr<SystemState<D>> systemPtr_in)
                 : settingsPtr(settingsPtr_in),
                   systemPtr(systemPtr_in){};

template <unsigned int D>
void Evolver<D>::allocate_cache()
{
    n_particles = systemPtr->n_particles;

    //Allocate on device memory for the system state before
    //evolution begins.
    evolver_cache = Cabana::AoSoA<EvolverCache, DeviceType, VECTOR_LENGTH>("cache", n_particles);
}

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

//==============================================================================
// this routine is used to initialize quantities prior to RK evolution
template <unsigned int D>
void Evolver<D>::set_current_timestep_quantities()
{

  CREATE_VIEW(device_, systemPtr->cabana_particles);

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
  Cabana::simd_parallel_for(*(systemPtr->simd_policy), fill_cache, "fill_cache");
}

//==========================================================================
template <unsigned int D>
void Evolver<D>:: advance_timestep_rk2( double dt,
                                        std::function<void(void)> time_derivatives_functional )
{
      double t0      = systemPtr->t;
      //double E0      = systemPtr->Ez;
      Kokkos::View<int, MemorySpace> E0("E0");
      Kokkos::View<int, MemorySpace> E("E");
      Kokkos::deep_copy(E0, systemPtr->Ez);
      Kokkos::deep_copy(E, systemPtr->Ez);
      // initialize quantities at current time step
      set_current_timestep_quantities();

      ////////////////////////////////////////////
      //    first step
      ////////////////////////////////////////////
      formatted_output::report("RK(n=2) evolution, step 1");

      // compute derivatives
      time_derivatives_functional();

  	  //Create views for the device
      CREATE_VIEW(device_, systemPtr->cabana_particles);
      auto slice_shv0 = Cabana::slice<evolver_cache_info::viscous_shear>(evolver_cache);
      auto slice_u0 = Cabana::slice<evolver_cache_info::four_velocity>(evolver_cache);
      auto slice_r0 = Cabana::slice<evolver_cache_info::position>(evolver_cache);
      auto slice_s0 = Cabana::slice<evolver_cache_info::specific_entropy>(evolver_cache);
      auto slice_Bulk0 = Cabana::slice<evolver_cache_info::Bulk_pressure>(evolver_cache);
      //auto E0 = Cabana::slice<evolver_cache_info::E0>(evolver_cache);

      auto update_rk2_step1 = KOKKOS_CLASS_LAMBDA(const int is, const int ia){
        
        //Cache previous step
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
        double specific_s = specific_s0 + .5*dt*d_dt_specific_s;
        double Bulk = Bulk0 + .5*dt*dBulk_dt;
        milne::Vector<double,D> r = r0 + .5*dt*v;
        milne::Vector<double,D> u = u0 + .5*dt*du_dt;
        milne::Matrix<double,D,D> shv;
        for(int idir=1; idir<D+1; ++idir)
        for(int jdir=1; jdir<D+1; ++jdir)
          shv(idir,jdir) = shv0(idir,jdir) + .5*dt*dshv_dt(idir-1, jdir-1);

        //Regulate negative entropy in edges
        if ( specific_s < 0.0 && specific_s0 < 5.E-1 ) specific_s = .5*(specific_s0+1);
        
        // regulate updated results if necessary
        if ( REGULATE_LARGE_S && specific_s > 10.0*specific_s0 ) specific_s = 2.0*specific_s0;

        //Update in memory
        device_specific_density.access(is, ia, densities_info::s) = specific_s;
        device_hydro_scalar.access(is, ia, hydro_info::Bulk) = Bulk;
        for (int idir=0; idir<D; ++idir){
          device_position.access(is, ia, idir) = r(idir);
          device_hydro_vector.access(is, ia, hydro_info::u, idir) = u(idir);
          for (int jdir=1; jdir<D; ++jdir)
            device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir) = shv(idir,jdir);
        }
      };
      Cabana::simd_parallel_for(*(systemPtr->simd_policy), update_rk2_step1, "update_rk2_step1");
      Kokkos::fence();
      
      systemPtr->t  = t0 + 0.5*dt;
      //Update hydro time on device
      double t = systemPtr->t;
      Cabana::simd_parallel_for(*(systemPtr->simd_policy), KOKKOS_CLASS_LAMBDA(const int is, const int ia)
      {
        device_hydro_scalar.access(is, ia, ccake::hydro_info::t) = t;
      }, "update_hydro_time_step1");
      Kokkos::fence();

      ////////////////////////////////////////////
      //    second step
      ////////////////////////////////////////////
      formatted_output::report("RK(n=2) evolution, step 2");

      // compute derivatives
      time_derivatives_functional();

      auto update_rk2_step2 = KOKKOS_CLASS_LAMBDA(const int is, const int ia){
        
        //Cache previous step
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
        double specific_s = specific_s0 + dt*d_dt_specific_s;
        double Bulk = Bulk0 + dt*dBulk_dt;
        milne::Vector<double,D> r = r0 + dt*v;
        milne::Vector<double,D> u = u0 + dt*du_dt;
        milne::Matrix<double,D,D> shv;
        for(int idir=1; idir<D+1; ++idir)
        for(int jdir=1; jdir<D+1; ++jdir)
          shv(idir,jdir) = shv0(idir,jdir) + dt*dshv_dt(idir-1, jdir-1);

        //Regulate negative entropy in edges
        if ( specific_s < 0.0 && specific_s0 < 5.E-1 ) specific_s = .5*(specific_s0+1);
        
        // regulate updated results if necessary
        if ( REGULATE_LARGE_S && specific_s > 10.0*specific_s0 ) specific_s = 2.0*specific_s0;

        //Update in memory
        device_specific_density.access(is, ia, densities_info::s) = specific_s;
        device_hydro_scalar.access(is, ia, hydro_info::Bulk) = Bulk;
        for (int idir=0; idir<D; ++idir){
          device_position.access(is, ia, idir) = r(idir);
          device_hydro_vector.access(is, ia, hydro_info::u, idir) = u(idir);
          for (int jdir=1; jdir<D; ++jdir)
            device_hydro_spacetime_matrix.access(is, ia, hydro_info::shv, idir, jdir) = shv(idir,jdir);
        }
      };
      Cabana::simd_parallel_for(*(systemPtr->simd_policy), update_rk2_step2, "update_rk2_step2");
      Kokkos::fence();
      systemPtr->t  = t0 + dt;

      //Update hydro time on device
      t = systemPtr->t;
      Cabana::simd_parallel_for(*(systemPtr->simd_policy), KOKKOS_CLASS_LAMBDA(const int is, const int ia)
      {
        device_hydro_scalar.access(is, ia,ccake::hydro_info::t) = t;
      }, "update_hydro_time_step2");
      Kokkos::fence();


      return;
}