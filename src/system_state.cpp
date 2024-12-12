#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>

#include "system_state.h"
#include "utilities.h"

using namespace constants;
using namespace ccake;

//Template instantiations
template class SystemState<1>;
template class SystemState<2>;
template class SystemState<3>;


/// @brief Initializes the system state.
/// @details Time is set to the initial state and smoothing lenght stored from
/// the settings.
/// @tparam D
template<unsigned int D>
void SystemState<D>::initialize()  // formerly called "manualenter"
{
  formatted_output::report("Initializing system");
  number_of_elapsed_timesteps = 0;
  t = settingsPtr->t0;
  return;
}

/// @brief Allocates memory in device for the particles and copies the data.
/// @see copy_host_to_device
template<unsigned int D>
void SystemState<D>::allocate_cabana_particles(){
    formatted_output::detail("Initializing device memory");
    ///Allocate memory for the particles
    cabana_particles = Cabana::AoSoA<CabanaParticle, DeviceType, VECTOR_LENGTH>("particles", n_particles);
    copy_host_to_device();
}

/// @brief Copies all the data from the array of particles to the AoSoA in device memory.
/// @tparam D The dimensionality of the simulation
/// @details This function creates an auxiliary AoSoA in host memory space and
/// fills it with the data from the particles array. Then it copies the data to
/// the device memory using a deep copy. Beware that this function is expensive
/// computationally and should be used only if strictly necessary.
template<unsigned int D>
void SystemState<D>::copy_host_to_device(){

  using SerialHost = Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>;
  //Auxiliary AoSoA for copying particles from/to host
  Cabana::AoSoA<ccake::CabanaParticle, SerialHost, 8> particles_h("particles_h",n_particles);
  #ifdef DEBUG
  formatted_output::detail("Copying data from host to device");
  #endif
  //Create vies to particles
  CREATE_VIEW(host_, particles_h);
  //Fill host arrays
  for (int iparticle=0; iparticle < n_particles; ++iparticle){
    //Copy hydro space matrix
    for (int i=0; i<D; ++i)
    for (int j=0; j<D; ++j){
        host_hydro_space_matrix(iparticle, ccake::hydro_info::Imat, i, j) = particles[iparticle].hydro.Imat(i, j);
        host_hydro_space_matrix(iparticle, ccake::hydro_info::gradV, i, j) = particles[iparticle].hydro.gradV(i, j);
        host_hydro_space_matrix(iparticle, ccake::hydro_info::M_u, i, j) = particles[iparticle].hydro.M_u(i, j);
        host_hydro_space_matrix(iparticle, ccake::hydro_info::R_u, i, j) = particles[iparticle].hydro.R_u(i, j);
        host_hydro_space_matrix(iparticle, ccake::hydro_info::R_0i_shear, i, j) = particles[iparticle].hydro.R_0i_shear(i, j);
        host_hydro_space_matrix(iparticle, ccake::hydro_info::M_0i_shear, i, j) = particles[iparticle].hydro.M_0i_shear(i, j);
      #ifdef DEBUG_SLOW
      //Couting to check values
      if (iparticle > 998 and iparticle < 1002) {
        std::cout << "gradV: " << particles[iparticle].hydro.gradV(i,j)  << std::endl;
      }
      #endif
    }
    host_hydro_scalar(iparticle, ccake::hydro_info::t) = particles[iparticle].hydro.t;
    host_hydro_scalar(iparticle, ccake::hydro_info::bulk) = particles[iparticle].hydro.bulk;
    host_hydro_scalar(iparticle, ccake::hydro_info::bigBulk) = particles[iparticle].hydro.bigBulk;
    host_hydro_scalar(iparticle, ccake::hydro_info::a) = particles[iparticle].hydro.a;
    host_hydro_scalar(iparticle, ccake::hydro_info::rho_B_ext) = particles[iparticle].hydro.rho_B_ext;
    host_hydro_scalar(iparticle, ccake::hydro_info::rho_S_ext) = particles[iparticle].hydro.rho_S_ext;
    host_hydro_scalar(iparticle, ccake::hydro_info::rho_Q_ext) = particles[iparticle].hydro.rho_Q_ext;
    host_hydro_scalar(iparticle, ccake::hydro_info::tau_Pi) = particles[iparticle].hydro.tau_Pi;
    host_hydro_scalar(iparticle, ccake::hydro_info::tau_pi) = particles[iparticle].hydro.tau_pi;
    host_hydro_scalar(iparticle, ccake::hydro_info::zeta_Pi) = particles[iparticle].hydro.zeta_Pi;
    host_hydro_scalar(iparticle, ccake::hydro_info::sigma_star) = particles[iparticle].hydro.sigma_star;
    host_hydro_scalar(iparticle, ccake::hydro_info::sigma) = particles[iparticle].hydro.sigma;
    host_hydro_scalar(iparticle, ccake::hydro_info::gamma) = particles[iparticle].hydro.gamma;
    host_hydro_scalar(iparticle, ccake::hydro_info::theta) = particles[iparticle].hydro.theta;
    host_hydro_scalar(iparticle, ccake::hydro_info::delta_PiPi) = particles[iparticle].hydro.delta_PiPi;
    host_hydro_scalar(iparticle, ccake::hydro_info::delta_pipi) = particles[iparticle].hydro.delta_pipi;
    host_hydro_scalar(iparticle, ccake::hydro_info::lambda_Pipi) = particles[iparticle].hydro.lambda_Pipi;
    host_hydro_scalar(iparticle, ccake::hydro_info::phi1) = particles[iparticle].hydro.phi1;
    host_hydro_scalar(iparticle, ccake::hydro_info::phi3) = particles[iparticle].hydro.phi3;
    host_hydro_scalar(iparticle, ccake::hydro_info::phi6) = particles[iparticle].hydro.phi6;
    host_hydro_scalar(iparticle, ccake::hydro_info::phi7) = particles[iparticle].hydro.phi7;
    host_hydro_scalar(iparticle, ccake::hydro_info::eta_pi) = particles[iparticle].hydro.eta_pi;
    host_hydro_scalar(iparticle, ccake::hydro_info::tau_pipi) = particles[iparticle].hydro.tau_pipi;
    host_hydro_scalar(iparticle, ccake::hydro_info::lambda_piPi) = particles[iparticle].hydro.lambda_piPi;
    host_hydro_scalar(iparticle, ccake::hydro_info::varsigma) = particles[iparticle].hydro.varsigma;
    host_hydro_scalar(iparticle, ccake::hydro_info::dbigBulk_dt) = particles[iparticle].hydro.dbigBulk_dt;
    host_hydro_scalar(iparticle, ccake::hydro_info::F_big_bulk) = particles[iparticle].hydro.F_big_bulk;
    host_hydro_scalar(iparticle, ccake::hydro_info::F_shv_nabla_u) = particles[iparticle].hydro.F_shv_nabla_u;
    host_hydro_scalar(iparticle, ccake::hydro_info::F_big_entropy) = particles[iparticle].hydro.F_big_entropy;
    host_hydro_scalar(iparticle, ccake::hydro_info::shv_nabla_u) = particles[iparticle].hydro.shv_nabla_u;


    

    // Copy hydro vector quantities
    for (int i = 0; i < D; ++i) {
        host_hydro_vector(iparticle, ccake::hydro_info::v, i) = particles[iparticle].hydro.v(i);
        host_hydro_vector(iparticle, ccake::hydro_info::u, i) = particles[iparticle].hydro.u(i);
        host_hydro_vector(iparticle, ccake::hydro_info::gradP, i) = particles[iparticle].hydro.gradP(i);
        host_hydro_vector(iparticle, ccake::hydro_info::gradE, i) = particles[iparticle].hydro.gradE(i);
        host_hydro_vector(iparticle, ccake::hydro_info::gradBulk, i) = particles[iparticle].hydro.gradBulk(i);
        host_hydro_vector(iparticle, ccake::hydro_info::divshear, i) = particles[iparticle].hydro.divshear(i);
        host_hydro_vector(iparticle, ccake::hydro_info::gradshear, i) = particles[iparticle].hydro.gradshear(i);
        host_hydro_vector(iparticle, ccake::hydro_info::du_dt, i) = particles[iparticle].hydro.du_dt(i);

        host_hydro_vector(iparticle, ccake::hydro_info::M_big_bulk, i) = particles[iparticle].hydro.M_big_bulk(i);
        host_hydro_vector(iparticle, ccake::hydro_info::M_shv_nabla_u, i) = particles[iparticle].hydro.M_shv_nabla_u(i);
        host_hydro_vector(iparticle, ccake::hydro_info::M_big_entropy, i) = particles[iparticle].hydro.M_big_entropy(i);
        host_hydro_vector(iparticle, ccake::hydro_info::F_u, i) = particles[iparticle].hydro.F_u(i);
       
        host_hydro_vector(iparticle, ccake::hydro_info::F_0i_shear, i) = particles[iparticle].hydro.F_0i_shear(i);
        host_hydro_vector(iparticle, ccake::hydro_info::j_ext, i) = particles[iparticle].hydro.j_ext(i);
    
        host_position(iparticle, i) = particles[iparticle].r(i);
    }
    //charges
    for(int i=0; i<3; ++i){
      host_hydro_vector(iparticle, ccake::hydro_info::F_big_N, i) = particles[iparticle].hydro.F_big_N(i);
      host_hydro_vector(iparticle, ccake::hydro_info::R_big_entropy, i) = particles[iparticle].hydro.R_big_entropy(i);
      host_hydro_vector(iparticle, ccake::hydro_info::R_big_bulk, i) = particles[iparticle].hydro.R_big_bulk(i);
      for(int j=0; j<3; ++j){
        host_hydro_space_matrix(iparticle, ccake::hydro_info::R_big_N, i, j) = particles[iparticle].hydro.R_big_N(i, j);
      }
      for(int idir=0; idir<D; ++idir){
        host_hydro_space_matrix(iparticle, ccake::hydro_info::M_big_N, idir, i) = particles[iparticle].hydro.M_big_N(idir, i);
      }
    }

    for (int i=D; i<3; i++) host_position(iparticle, i) = 0;

    host_thermo(iparticle, ccake::thermo_info::T) = particles[iparticle].thermo.T;
    host_thermo(iparticle, ccake::thermo_info::muB) = particles[iparticle].thermo.muB;
    host_thermo(iparticle, ccake::thermo_info::muS) = particles[iparticle].thermo.muS;
    host_thermo(iparticle, ccake::thermo_info::muQ) = particles[iparticle].thermo.muQ;
    host_thermo(iparticle, ccake::thermo_info::e) = particles[iparticle].thermo.e;
    host_thermo(iparticle, ccake::thermo_info::s) = particles[iparticle].thermo.s;
    host_thermo(iparticle, ccake::thermo_info::rhoB) = particles[iparticle].thermo.rhoB;
    host_thermo(iparticle, ccake::thermo_info::rhoS) = particles[iparticle].thermo.rhoS;
    host_thermo(iparticle, ccake::thermo_info::rhoQ) = particles[iparticle].thermo.rhoQ;
    host_thermo(iparticle, ccake::thermo_info::p) = particles[iparticle].thermo.p;
    host_thermo(iparticle, ccake::thermo_info::cs2) = particles[iparticle].thermo.cs2;
    host_thermo(iparticle, ccake::thermo_info::w) = particles[iparticle].thermo.w;
    host_thermo(iparticle, ccake::thermo_info::A) = particles[iparticle].thermo.A;
    host_thermo(iparticle, ccake::thermo_info::dwds) = particles[iparticle].thermo.dwds;
    host_thermo(iparticle, ccake::thermo_info::dwdB) = particles[iparticle].thermo.dwdB;
    host_thermo(iparticle, ccake::thermo_info::dwdS) = particles[iparticle].thermo.dwdS;
    host_thermo(iparticle, ccake::thermo_info::dwdQ) = particles[iparticle].thermo.dwdQ;
    #ifdef DEBUG_SLOW
    //Couting to check values
    if (iparticle > 998 and iparticle < 1002) {
      std::cout << " t " << particles[iparticle].hydro.t << " gamma " << particles[iparticle].hydro.gamma << \
              " sigma_star " << particles[iparticle].hydro.sigma_star << " T " << particles[iparticle].thermo.T \
	     << " thermo.e " << particles[iparticle].thermo.e << " thermo.s " << particles[iparticle].thermo.s << std::endl;
    }
    #endif

  // Copy hydro spacetime matrix
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            host_hydro_spacetime_matrix(iparticle, ccake::hydro_info::shv, i, j) = particles[iparticle].hydro.shv(i, j);

    // Copy hydro shear auxiliary vectors and matrices
    for (int i = 0; i < 2; ++i) {
        for (int j = i; j < 3; ++j) {
            host_hydro_shear_aux_vector(iparticle, ccake::hydro_info::F_big_shear, i, j) = particles[iparticle].hydro.F_big_shear(i, j);
            host_hydro_shear_aux_vector(iparticle, ccake::hydro_info::d_bigshv_dt, i, j) = particles[iparticle].hydro.d_bigshv_dt(i, j);
            host_hydro_shear_aux_vector(iparticle, ccake::hydro_info::bigshv, i, j) = particles[iparticle].hydro.bigshv(i, j);
            for (int k = 0; k < D; ++k) {
                int linear_index = j * D + i;
                host_hydro_shear_aux_matrix(iparticle, ccake::hydro_info::M_big_shear, i, linear_index) = particles[iparticle].hydro.M_big_shear(i, linear_index);
            }
            for (int k = 0; k < 3; ++k) {
                int linear_index = j * 3 + i;
                host_hydro_shear_aux_matrix(iparticle, ccake::hydro_info::R_big_shear, i, linear_index) = particles[iparticle].hydro.R_big_shear(i, linear_index);
            }
        }
    }

    host_input(iparticle, ccake::densities_info::e) = particles[iparticle].input.e;
    host_input(iparticle, ccake::densities_info::s) = particles[iparticle].input.s;
    host_input(iparticle, ccake::densities_info::rhoB) = particles[iparticle].input.rhoB;
    host_input(iparticle, ccake::densities_info::rhoS) = particles[iparticle].input.rhoS;
    host_input(iparticle, ccake::densities_info::rhoQ) = particles[iparticle].input.rhoQ;

    host_smoothed(iparticle, ccake::densities_info::e) = particles[iparticle].smoothed.e;
    host_smoothed(iparticle, ccake::densities_info::s) = particles[iparticle].smoothed.s;
    host_smoothed(iparticle, ccake::densities_info::rhoB) = particles[iparticle].smoothed.rhoB;
    host_smoothed(iparticle, ccake::densities_info::rhoS) = particles[iparticle].smoothed.rhoS;
    host_smoothed(iparticle, ccake::densities_info::rhoQ) = particles[iparticle].smoothed.rhoQ;
    #ifdef DEBUG_SLOW
    //Couting to check values
    if (iparticle > 998 and iparticle < 1002) {
      std::cout << " smoothed.e " << particles[iparticle].smoothed.e << " smoothed.s " << particles[iparticle].smoothed.s << std::endl;
			      }
    #endif

    host_specific_density(iparticle, ccake::densities_info::e) = particles[iparticle].specific.e;
    host_specific_density(iparticle, ccake::densities_info::s) = particles[iparticle].specific.s;
    host_specific_density(iparticle, ccake::densities_info::rhoB) = particles[iparticle].specific.rhoB;
    host_specific_density(iparticle, ccake::densities_info::rhoS) = particles[iparticle].specific.rhoS;
    host_specific_density(iparticle, ccake::densities_info::rhoQ) = particles[iparticle].specific.rhoQ;

    host_d_dt_spec(iparticle, ccake::densities_info::e) = particles[iparticle].d_dt_spec.e;
    host_d_dt_spec(iparticle, ccake::densities_info::s) = particles[iparticle].d_dt_spec.s;
    host_d_dt_spec(iparticle, ccake::densities_info::rhoB) = particles[iparticle].d_dt_spec.rhoB;
    host_d_dt_spec(iparticle, ccake::densities_info::rhoS) = particles[iparticle].d_dt_spec.rhoS;
    host_d_dt_spec(iparticle, ccake::densities_info::rhoQ) = particles[iparticle].d_dt_spec.rhoQ;

    host_norm_spec(iparticle, ccake::densities_info::e) = particles[iparticle].norm_spec.e;
    host_norm_spec(iparticle, ccake::densities_info::s) = particles[iparticle].norm_spec.s;
    host_norm_spec(iparticle, ccake::densities_info::rhoB) = particles[iparticle].norm_spec.rhoB;
    host_norm_spec(iparticle, ccake::densities_info::rhoS) = particles[iparticle].norm_spec.rhoS;
    host_norm_spec(iparticle, ccake::densities_info::rhoQ) = particles[iparticle].norm_spec.rhoQ;

    host_efcheck(iparticle) = particles[iparticle].efcheck;
    host_contribution_to_total_E(iparticle)  = particles[iparticle].contribution_to_total_E;
    host_contribution_to_total_dEz(iparticle)  = particles[iparticle].contribution_to_total_dEz;
    host_contribution_to_total_Ez(iparticle)  = particles[iparticle].contribution_to_total_Ez;
    host_id(iparticle)  = particles[iparticle].ID;
    host_btrack(iparticle)  = particles[iparticle].btrack;
    host_freeze(iparticle)  = particles[iparticle].Freeze;
  }
  Cabana::deep_copy( cabana_particles, particles_h );
  Kokkos::fence();

}

/// @brief Takes the data out of the AoSoA in device memory back to particle
/// data structure
/// @details This function creates a temporary AoSoA in host memory space and
/// mirrors the device memory. Just then it populates the particles arrays with
/// the data from the host AoSoA. Beware that this process is expensive and
/// should be used only if strictly necessary.
template<unsigned int D>
void SystemState<D>::copy_device_to_host(){

  #ifdef DEBUG
  formatted_output::detail("Copying data from device to host");
  #endif
  //Auxiliary AoSoA for copying particles from/to host
  auto particles_h = Cabana::create_mirror_view_and_copy(Kokkos::HostSpace(), cabana_particles);
  Kokkos::fence();
  CREATE_VIEW(host_, particles_h);
  for (int iparticle=0; iparticle < n_particles; ++iparticle){
    int id = host_id(iparticle);
    //Copy hydro space matrix
        for (int i=0; i<D; ++i)
    for (int j=0; j<D; ++j){
      particles[id].hydro.Imat(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::Imat, i, j);
      particles[id].hydro.gradV(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::gradV, i, j);
      particles[id].hydro.M_u(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::M_u, i, j);
      particles[id].hydro.R_u(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::R_u, i, j);
      particles[id].hydro.R_0i_shear(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::R_0i_shear, i, j);
      particles[id].hydro.M_0i_shear(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::M_0i_shear, i, j);
    }
    particles[id].hydro.t = host_hydro_scalar(iparticle, ccake::hydro_info::t);
    particles[id].hydro.bulk = host_hydro_scalar(iparticle, ccake::hydro_info::bulk);
    particles[id].hydro.bigBulk = host_hydro_scalar(iparticle, ccake::hydro_info::bigBulk);
    particles[id].hydro.a = host_hydro_scalar(iparticle, ccake::hydro_info::a);
    particles[id].hydro.rho_B_ext = host_hydro_scalar(iparticle, ccake::hydro_info::rho_B_ext);
    particles[id].hydro.rho_S_ext = host_hydro_scalar(iparticle, ccake::hydro_info::rho_S_ext);
    particles[id].hydro.rho_Q_ext = host_hydro_scalar(iparticle, ccake::hydro_info::rho_Q_ext);
    particles[id].hydro.tau_Pi = host_hydro_scalar(iparticle, ccake::hydro_info::tau_Pi);
    particles[id].hydro.tau_pi = host_hydro_scalar(iparticle, ccake::hydro_info::tau_pi);
    particles[id].hydro.zeta_Pi = host_hydro_scalar(iparticle, ccake::hydro_info::zeta_Pi);
    particles[id].hydro.sigma_star = host_hydro_scalar(iparticle, ccake::hydro_info::sigma_star);
    particles[id].hydro.sigma = host_hydro_scalar(iparticle, ccake::hydro_info::sigma);
    particles[id].hydro.gamma = host_hydro_scalar(iparticle, ccake::hydro_info::gamma);
    particles[id].hydro.theta = host_hydro_scalar(iparticle, ccake::hydro_info::theta);
    particles[id].hydro.delta_PiPi = host_hydro_scalar(iparticle, ccake::hydro_info::delta_PiPi);
    particles[id].hydro.delta_pipi = host_hydro_scalar(iparticle, ccake::hydro_info::delta_pipi);
    particles[id].hydro.lambda_Pipi = host_hydro_scalar(iparticle, ccake::hydro_info::lambda_Pipi);
    particles[id].hydro.phi1 = host_hydro_scalar(iparticle, ccake::hydro_info::phi1);
    particles[id].hydro.phi3 = host_hydro_scalar(iparticle, ccake::hydro_info::phi3);
    particles[id].hydro.phi6 = host_hydro_scalar(iparticle, ccake::hydro_info::phi6);
    particles[id].hydro.phi7 = host_hydro_scalar(iparticle, ccake::hydro_info::phi7);
    particles[id].hydro.eta_pi = host_hydro_scalar(iparticle, ccake::hydro_info::eta_pi);
    particles[id].hydro.tau_pipi = host_hydro_scalar(iparticle, ccake::hydro_info::tau_pipi);
    particles[id].hydro.lambda_piPi = host_hydro_scalar(iparticle, ccake::hydro_info::lambda_piPi);
    particles[id].hydro.varsigma = host_hydro_scalar(iparticle, ccake::hydro_info::varsigma);
    particles[id].hydro.dbigBulk_dt = host_hydro_scalar(iparticle, ccake::hydro_info::dbigBulk_dt);
    particles[id].hydro.F_big_bulk = host_hydro_scalar(iparticle, ccake::hydro_info::F_big_bulk);
    particles[id].hydro.F_shv_nabla_u = host_hydro_scalar(iparticle, ccake::hydro_info::F_shv_nabla_u);
    particles[id].hydro.F_big_entropy = host_hydro_scalar(iparticle, ccake::hydro_info::F_big_entropy);
    particles[id].hydro.shv_nabla_u = host_hydro_scalar(iparticle, ccake::hydro_info::shv_nabla_u);

    for (int i=0; i<D; i++){
      particles[id].hydro.v(i) = host_hydro_vector(iparticle, ccake::hydro_info::v,i);
      particles[id].hydro.u(i) = host_hydro_vector(iparticle, ccake::hydro_info::u,i);
      particles[id].hydro.gradP(i) = host_hydro_vector(iparticle, ccake::hydro_info::gradP,i);
      particles[id].hydro.gradE(i) = host_hydro_vector(iparticle, ccake::hydro_info::gradE,i);
      particles[id].hydro.gradBulk(i) = host_hydro_vector(iparticle, ccake::hydro_info::gradBulk,i);
      particles[id].hydro.divshear(i) = host_hydro_vector(iparticle, ccake::hydro_info::divshear,i);
      particles[id].hydro.gradshear(i) = host_hydro_vector(iparticle, ccake::hydro_info::gradshear,i);
      particles[id].hydro.du_dt(i) = host_hydro_vector(iparticle, ccake::hydro_info::du_dt,i);
      particles[id].hydro.M_big_bulk(i) = host_hydro_vector(iparticle, ccake::hydro_info::M_big_bulk,i);
      particles[id].hydro.M_shv_nabla_u(i) = host_hydro_vector(iparticle, ccake::hydro_info::M_shv_nabla_u,i);
      particles[id].hydro.M_big_entropy(i) = host_hydro_vector(iparticle, ccake::hydro_info::M_big_entropy,i);
      particles[id].hydro.F_u(i) = host_hydro_vector(iparticle, ccake::hydro_info::F_u,i);
      particles[id].hydro.F_0i_shear(i) = host_hydro_vector(iparticle, ccake::hydro_info::F_0i_shear,i);
      particles[id].hydro.j_ext(i) = host_hydro_vector(iparticle, ccake::hydro_info::j_ext,i);
      particles[id].r(i) = host_position(iparticle, i);
    }
    //charges
    for(int i=0; i<3; ++i){
      particles[id].hydro.F_big_N(i) = host_hydro_vector(iparticle, ccake::hydro_info::F_big_N,i);
      particles[id].hydro.R_big_entropy(i) = host_hydro_vector(iparticle, ccake::hydro_info::R_big_entropy,i);
      particles[id].hydro.R_big_bulk(i) = host_hydro_vector(iparticle, ccake::hydro_info::R_big_bulk,i);
      for(int j=0; j<3; ++j){
        particles[id].hydro.R_big_N(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::R_big_N,i,j);
      }
      for(int idir=0; idir<D; idir++){
        particles[id].hydro.M_big_N(idir,i) = host_hydro_space_matrix(iparticle, ccake::hydro_info::M_big_N,idir,i);
      }
    }

    particles[id].thermo.T = host_thermo(iparticle, ccake::thermo_info::T);
    particles[id].thermo.muB = host_thermo(iparticle, ccake::thermo_info::muB);
    particles[id].thermo.muS = host_thermo(iparticle, ccake::thermo_info::muS);
    particles[id].thermo.muQ = host_thermo(iparticle, ccake::thermo_info::muQ);
    particles[id].thermo.e = host_thermo(iparticle, ccake::thermo_info::e);
    particles[id].thermo.s = host_thermo(iparticle, ccake::thermo_info::s);
    particles[id].thermo.rhoB = host_thermo(iparticle, ccake::thermo_info::rhoB);
    particles[id].thermo.rhoS = host_thermo(iparticle, ccake::thermo_info::rhoS);
    particles[id].thermo.rhoQ = host_thermo(iparticle, ccake::thermo_info::rhoQ);
    particles[id].thermo.p = host_thermo(iparticle, ccake::thermo_info::p);
    particles[id].thermo.cs2 = host_thermo(iparticle, ccake::thermo_info::cs2);
    particles[id].thermo.w = host_thermo(iparticle, ccake::thermo_info::w);
    particles[id].thermo.A = host_thermo(iparticle, ccake::thermo_info::A);
    particles[id].thermo.dwds = host_thermo(iparticle, ccake::thermo_info::dwds);
    particles[id].thermo.dwdB = host_thermo(iparticle, ccake::thermo_info::dwdB);
    particles[id].thermo.dwdS = host_thermo(iparticle, ccake::thermo_info::dwdS);
    particles[id].thermo.dwdQ = host_thermo(iparticle, ccake::thermo_info::dwdQ);

    for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      particles[id].hydro.shv(i,j) = host_hydro_spacetime_matrix(iparticle, ccake::hydro_info::shv, i, j);
    for (int i=0; i<2; i++){
      for (int j=i; j<3; j++){
        particles[id].hydro.F_big_shear(i,j) = host_hydro_shear_aux_vector(iparticle, ccake::hydro_info::F_big_shear, i, j);
        particles[id].hydro.d_bigshv_dt(i,j) = host_hydro_shear_aux_vector(iparticle, ccake::hydro_info::d_bigshv_dt, i, j);
        particles[id].hydro.bigshv(i,j) = host_hydro_shear_aux_vector(iparticle, ccake::hydro_info::bigshv, i, j);
        for (int k=0; k<D; k++){
          int linear_index = j * D + k;
          particles[id].hydro.M_big_shear(i,linear_index) = host_hydro_shear_aux_matrix(iparticle, ccake::hydro_info::M_big_shear, i, linear_index);
        }
        for (int k=0; k<3; k++){
          int linear_index = j * 3 + k;
          particles[id].hydro.R_big_shear(i,linear_index) = host_hydro_shear_aux_matrix(iparticle, ccake::hydro_info::R_big_shear, i, linear_index);
        }
      }
    }
      

    particles[id].input.e = host_input(iparticle, ccake::densities_info::e);
    particles[id].input.s = host_input(iparticle, ccake::densities_info::s);
    particles[id].input.rhoB = host_input(iparticle, ccake::densities_info::rhoB);
    particles[id].input.rhoS = host_input(iparticle, ccake::densities_info::rhoS);
    particles[id].input.rhoQ = host_input(iparticle, ccake::densities_info::rhoQ);

    particles[id].smoothed.e = host_smoothed(iparticle, ccake::densities_info::e);
    particles[id].smoothed.s = host_smoothed(iparticle, ccake::densities_info::s);
    particles[id].smoothed.rhoB = host_smoothed(iparticle, ccake::densities_info::rhoB);
    particles[id].smoothed.rhoS = host_smoothed(iparticle, ccake::densities_info::rhoS);
    particles[id].smoothed.rhoQ = host_smoothed(iparticle, ccake::densities_info::rhoQ);

    particles[id].specific.e = host_specific_density(iparticle, ccake::densities_info::e);
    particles[id].specific.s = host_specific_density(iparticle, ccake::densities_info::s);
    particles[id].specific.rhoB = host_specific_density(iparticle, ccake::densities_info::rhoB);
    particles[id].specific.rhoS = host_specific_density(iparticle, ccake::densities_info::rhoS);
    particles[id].specific.rhoQ = host_specific_density(iparticle, ccake::densities_info::rhoQ);

    particles[id].d_dt_spec.e = host_d_dt_spec(iparticle, ccake::densities_info::e);
    particles[id].d_dt_spec.s = host_d_dt_spec(iparticle, ccake::densities_info::s);
    particles[id].d_dt_spec.rhoB = host_d_dt_spec(iparticle, ccake::densities_info::rhoB);
    particles[id].d_dt_spec.rhoS = host_d_dt_spec(iparticle, ccake::densities_info::rhoS);
    particles[id].d_dt_spec.rhoQ = host_d_dt_spec(iparticle, ccake::densities_info::rhoQ);

    particles[id].norm_spec.e = host_norm_spec(iparticle, ccake::densities_info::e);
    particles[id].norm_spec.s = host_norm_spec(iparticle, ccake::densities_info::s);
    particles[id].norm_spec.rhoB = host_norm_spec(iparticle, ccake::densities_info::rhoB);
    particles[id].norm_spec.rhoS = host_norm_spec(iparticle, ccake::densities_info::rhoS);
    particles[id].norm_spec.rhoQ = host_norm_spec(iparticle, ccake::densities_info::rhoQ);

    particles[id].efcheck = host_efcheck(iparticle);
    particles[id].contribution_to_total_E = host_contribution_to_total_E(iparticle);
    particles[id].contribution_to_total_dEz = host_contribution_to_total_dEz(iparticle);
    particles[id].contribution_to_total_Ez = host_contribution_to_total_Ez(iparticle);
    particles[id].ID = host_id(iparticle);
    particles[id].btrack = host_btrack(iparticle);
    particles[id].Freeze = host_freeze(iparticle);
  }

}

/// @brief Reset the neighbour list
/// @details This function should be called after the particles have been moved,
/// since neighbour list is based on the position of the particles. Computation
/// of nearest neighbor is done using an auxiliary grid to help find the
/// neighbours. The grid is defined as twice as the initial domain size. If the
/// particles are moved outside this domain, the grid should be updated or
/// a segmentation fault may occur.
///
/// Beware that finding nearest neighbor is an expensive operation and should be
/// be minimized. This makes it also a good candidate for optimization.
/// @todo We need to implement a way to update the grid when the particles are
/// moved outside the initial domain. Notice that the current setup incurs in
/// many empty cells which could cause a performance hit.
/// @todo To get the number of nearest neighbors, we need to use UVM memory.
/// This means an implicit copy from device memory to host memory is happening.
/// I don't know how bad this is affecting performance and should be further
/// investigated. Ideally, I would be happier if we could just get completelly
/// rid of UVM memory, since implicit copies could be introduced in other parts
/// of the code without us noticing.
/// @tparam D The dimensionality of the simulation
template<unsigned int D>
void SystemState<D>::reset_neighbour_list(){
  //Outfile for cabana particle neighbors
  /*ofstream outfile;
  outfile.open("neighbors.dat")*/

  double min_pos[3], max_pos[3];
  switch (D)
  {
    case 1:
      min_pos[0] = settingsPtr->etamin;
      min_pos[1] = settingsPtr->xmin;
      min_pos[2] = settingsPtr->ymin;
      break;
    default:
      min_pos[0] = settingsPtr->xmin;
      min_pos[1] = settingsPtr->ymin;
      min_pos[2] = settingsPtr->etamin;
      break;
  }

  for(int idir=0; idir<3; ++idir)
    min_pos[idir] *= 2.; //Grid must be 100% bigger ///TODO: Allow this to be an optional input parameter

  //Cabana needs a 3D grid. We set the remaining dimensions to be a single cell
  double neighborhood_radius = 2*settingsPtr->hT;
  for(int idir=D; idir<3; ++idir) min_pos[idir] = -settingsPtr->hT;
  for(int idir=0; idir<3; ++idir)
    max_pos[idir] = -min_pos[idir];

  CREATE_VIEW(device_, cabana_particles);
  //Enabling change the order of the particles in the AoSoA. This may be a problem.
  //Cabana::permute( cell_list, cabana_particles ); 
  double cell_ratio = 1.; //neighbour to cell_space ratio
  neighbour_list = ListType (  device_position, 0, device_position.size(),
                                          neighborhood_radius, cell_ratio, min_pos, max_pos
                           );
  Kokkos::fence();

  //Update the number of neighbours
  ///TODO: Requires UVM, which is not good for performance
  ///Need a way to call Cabana::NeighborList::numNeighbor directly from GPU
  ///without seg fault
  for (int i=0; i<n_particles; ++i) {
    device_btrack(i) = device_btrack(i) == -1 ? -1 : Cabana::NeighborList<ListType>::numNeighbor( neighbour_list, i );
    /*if (Cabana::NeighborList<ListType>::numNeighbor( neighbour_list, i ) < 5) {
      std::cout << "The particle " << i << " was found with less than 5 neighbors" << std::endl;
      abort();
    }*/
  }
  //print_neighbors(950);
  //print_neighbors(975);
  //print_neighbors(1000);
  //print_neighbors(1025);
  #ifdef DEBUG_SLOW
  print_neighbors(0);
  #endif
}

/// @brief Entry point for computing neighbor list. Name kept for compatibility.
/// @details This function is a wrapper for the actual computation of the
/// neighbour list. It is called by the main loop and should be called after the
/// particles have been moved.
/// @see reset_neighbour_list
template<unsigned int D>
void SystemState<D>::initialize_linklist()
{
  formatted_output::report("Initializing linklist");

  reset_neighbour_list();
  return;
}

/// @brief Computes total entropy
/// @tparam D The dimensionality of the simulation
/// @param first_iteration true if it is the first iteration
template<unsigned int D>
void SystemState<D>::conservation_entropy(bool first_iteration)
{

  CREATE_VIEW(device_, cabana_particles);
  S = 0.0;
  auto get_total_entropy = KOKKOS_LAMBDA(const int i, double &local_S)
  {
    local_S += device_specific_density(i, ccake::densities_info::s)*device_norm_spec(i, ccake::densities_info::s);
  };
  Kokkos::parallel_reduce("loop_conservation_entropy",n_particles, get_total_entropy, S);
  if (first_iteration) S0 = S;
}

/// @brief Computes total Baryon, Strangeness and Electric charge
/// @tparam D The dimensionality of the simulation
/// @param first_iteration true if it is the first iteration
template<unsigned int D>
void SystemState<D>::conservation_BSQ(bool first_iteration)
{
  // reset
  Btotal = 0.0;
  Stotal = 0.0;
  Qtotal = 0.0;
  CREATE_VIEW(device_, cabana_particles);
  auto get_total_B = KOKKOS_LAMBDA(const int &i, double &Btotal)
  {
    Btotal += device_specific_density(i, ccake::densities_info::rhoB)*device_norm_spec(i, ccake::densities_info::rhoB);
  };
  auto get_total_S = KOKKOS_LAMBDA(const int &i, double &Stotal)
  {
    Stotal += device_specific_density(i, ccake::densities_info::rhoS)*device_norm_spec(i, ccake::densities_info::rhoS);
  };
  auto get_total_Q = KOKKOS_LAMBDA(const int &i, double &Qtotal)
  {
    Qtotal += device_specific_density(i, ccake::densities_info::rhoQ)*device_norm_spec(i, ccake::densities_info::rhoQ);
  };
  Kokkos::parallel_reduce("loop_conservation_B",n_particles, get_total_B, Kokkos::Sum<double>(Btotal));
  Kokkos::parallel_reduce("loop_conservation_S",n_particles, get_total_S, Kokkos::Sum<double>(Stotal));
  Kokkos::parallel_reduce("loop_conservation_Q",n_particles, get_total_Q, Kokkos::Sum<double>(Qtotal));

  // save initial totals
  if (first_iteration)
  {
    Btotal0 = Btotal;
    Stotal0 = Stotal;
    Qtotal0 = Qtotal;
  }
  return;
}




///////////////////////////////////////
//TODO: Parallelize with Kokkos::parallel_reduce


template<unsigned int D>
void SystemState<D>::conservation_energy(bool first_iteration)
{
  ///////////////////////////////////////////////
  // don't bother checking energy conservation on
  // intermediate RK steps
  // calculate total energy (T^{00})
  double t = t;
  CREATE_VIEW(device_, cabana_particles);
  E = 0.0;
  auto get_total_energy = KOKKOS_LAMBDA(const int i, double &local_E)
  {
    double w = device_thermo(i, ccake::thermo_info::w);
    double dE = device_contribution_to_total_E(i);
    double p = device_thermo(i, ccake::thermo_info::p);  
    double norm_spec = device_norm_spec(i, ccake::densities_info::s);
    double sigma_star = device_hydro_scalar(i, ccake::hydro_info::sigma_star);
    double bulk = device_hydro_scalar(i, ccake::hydro_info::bulk);
    double shv00 = device_hydro_spacetime_matrix(i, ccake::hydro_info::shv, 0,0);
    double gamma = device_hydro_scalar(i, ccake::hydro_info::gamma);
    double g2 = gamma*gamma;

    double C = w + bulk;
    
    local_E += (C * g2 - p - bulk + shv00) * norm_spec * t / sigma_star;
  };
  Kokkos::parallel_reduce("loop_conservation_energy",n_particles, get_total_energy, E);
  Kokkos::fence();
  
  if (first_iteration) E0 = E;
  //calculate Ez

  Ez = 0.0;
  auto get_total_Ez = KOKKOS_LAMBDA(const int i, double &local_Ez)
  {
    local_Ez += device_contribution_to_total_Ez(i);
  };
  Kokkos::parallel_reduce("loop_conservation_Ez",n_particles, get_total_Ez, Ez);
  Kokkos::fence();

  Etot  = E + Ez;
  Eloss = (E0-Etot)/E0*100;
}

///////////////////////////////////////
template<unsigned int D>
void SystemState<D>::compute_eccentricities()
{
  timesteps.push_back( t );

  compute_e_2_X();
  compute_e_2_P();
}

///////////////////////////////////////
//TODO: Use Cabana for paralelism
template<unsigned int D>
void SystemState<D>::compute_e_2_P()
{
  double e_2_P_c = 0.0, e_2_P_s = 0.0, normalization = 0.0;

  for ( auto & p : particles )
  {
    double ux        = p.hydro.u(0),
           uy        = p.hydro.u(1);
    double p_plus_Pi = p.p() + p.hydro.bulk;
    double e_p_Pi    = p.e() + p_plus_Pi;

    double Txx = e_p_Pi*ux*ux + p_plus_Pi + p.hydro.shv(1,1);
    double Txy = e_p_Pi*ux*uy             + p.hydro.shv(1,2);
    double Tyy = e_p_Pi*uy*uy + p_plus_Pi + p.hydro.shv(2,2);

    e_2_P_c += Txx-Tyy;
    e_2_P_s += 2.0*Txy;
    normalization += Txx+Tyy;
  }
  e_2_P.push_back( sqrt(e_2_P_c*e_2_P_c+e_2_P_s*e_2_P_s)/abs(normalization) );
}

///////////////////////////////////////
template<unsigned int D>
//TODO: Use Cabana for paralelism
void SystemState<D>::compute_e_2_X()
{
  double e_2_X_c = 0.0, e_2_X_s = 0.0, normalization = 0.0;
  for ( auto & p : particles )
  {
    double x       = p.r(0),
           y       = p.r(1);
    e_2_X_c       += p.hydro.gamma*p.e()*(x*x-y*y);
    e_2_X_s       += 2.0*p.hydro.gamma*p.e()*x*y;
    normalization += p.hydro.gamma*p.e()*(x*x+y*y);
  }
  e_2_X.push_back( sqrt(e_2_X_c*e_2_X_c+e_2_X_s*e_2_X_s)/abs(normalization) );
}
