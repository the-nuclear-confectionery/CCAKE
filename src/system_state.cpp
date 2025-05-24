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
    host_hydro_scalar(iparticle, ccake::hydro_info::causality) = particles[iparticle].hydro.causality;
    host_hydro_scalar(iparticle, ccake::hydro_info::bulk) = particles[iparticle].hydro.bulk;
    host_hydro_scalar(iparticle, ccake::hydro_info::extensive_bulk) = particles[iparticle].hydro.extensive_bulk;
    host_hydro_scalar(iparticle, ccake::hydro_info::a) = particles[iparticle].hydro.a;
    host_hydro_scalar(iparticle, ccake::hydro_info::rho_B_ext) = particles[iparticle].hydro.rho_B_ext;
    host_hydro_scalar(iparticle, ccake::hydro_info::rho_S_ext) = particles[iparticle].hydro.rho_S_ext;
    host_hydro_scalar(iparticle, ccake::hydro_info::rho_Q_ext) = particles[iparticle].hydro.rho_Q_ext;
    host_hydro_scalar(iparticle, ccake::hydro_info::tau_Pi) = particles[iparticle].hydro.tau_Pi;
    host_hydro_scalar(iparticle, ccake::hydro_info::tau_pi) = particles[iparticle].hydro.tau_pi;
    host_hydro_scalar(iparticle, ccake::hydro_info::zeta_Pi) = particles[iparticle].hydro.zeta_Pi;
    host_hydro_scalar(iparticle, ccake::hydro_info::sigma_lab) = particles[iparticle].hydro.sigma_lab;
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
    host_hydro_scalar(iparticle, ccake::hydro_info::tmunu_trace) = particles[iparticle].hydro.tmunu_trace;
    host_hydro_scalar(iparticle, ccake::hydro_info::d_extensive_bulk_dt) = particles[iparticle].hydro.d_extensive_bulk_dt;
    host_hydro_scalar(iparticle, ccake::hydro_info::F_extensive_bulk) = particles[iparticle].hydro.F_extensive_bulk;
    host_hydro_scalar(iparticle, ccake::hydro_info::F_shv_nabla_u) = particles[iparticle].hydro.F_shv_nabla_u;
    host_hydro_scalar(iparticle, ccake::hydro_info::F_extensive_entropy) = particles[iparticle].hydro.F_extensive_entropy;
    host_hydro_scalar(iparticle, ccake::hydro_info::shv_nabla_u) = particles[iparticle].hydro.shv_nabla_u;
    host_hydro_scalar(iparticle, ccake::hydro_info::j0_ext) = particles[iparticle].hydro.j0_ext;
    host_hydro_scalar(iparticle, ccake::hydro_info::shear_knudsen) = particles[iparticle].hydro.shear_knudsen;
    host_hydro_scalar(iparticle, ccake::hydro_info::bulk_knudsen) = particles[iparticle].hydro.bulk_knudsen;
    host_hydro_scalar(iparticle, ccake::hydro_info::inverse_reynolds_shear) = particles[iparticle].hydro.inverse_reynolds_shear;
    host_hydro_scalar(iparticle, ccake::hydro_info::inverse_reynolds_bulk) = particles[iparticle].hydro.inverse_reynolds_bulk;



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

        host_hydro_vector(iparticle, ccake::hydro_info::M_extensive_bulk, i) = particles[iparticle].hydro.M_extensive_bulk(i);
        host_hydro_vector(iparticle, ccake::hydro_info::M_shv_nabla_u, i) = particles[iparticle].hydro.M_shv_nabla_u(i);
        host_hydro_vector(iparticle, ccake::hydro_info::M_extensive_entropy, i) = particles[iparticle].hydro.M_extensive_entropy(i);
        host_hydro_vector(iparticle, ccake::hydro_info::F_u, i) = particles[iparticle].hydro.F_u(i);

        host_hydro_vector(iparticle, ccake::hydro_info::F_0i_shear, i) = particles[iparticle].hydro.F_0i_shear(i);
        host_hydro_vector(iparticle, ccake::hydro_info::j_ext, i) = particles[iparticle].hydro.j_ext(i);

        host_position(iparticle, i) = particles[iparticle].r(i);
    }
    //charges
    for(int i=0; i<3; ++i){
      host_hydro_vector(iparticle, ccake::hydro_info::F_extensive_N, i) = particles[iparticle].hydro.F_extensive_N(i);
      host_hydro_vector(iparticle, ccake::hydro_info::R_extensive_entropy, i) = particles[iparticle].hydro.R_extensive_entropy(i);
      host_hydro_vector(iparticle, ccake::hydro_info::R_extensive_bulk, i) = particles[iparticle].hydro.R_extensive_bulk(i);
      for(int j=0; j<3; ++j){
        host_hydro_space_matrix(iparticle, ccake::hydro_info::R_extensive_N, i, j) = particles[iparticle].hydro.R_extensive_N(i, j);
      }
      for(int idir=0; idir<D; ++idir){
        host_hydro_space_matrix(iparticle, ccake::hydro_info::M_extensive_N, idir, i) = particles[iparticle].hydro.M_extensive_N(idir, i);
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
    host_thermo(iparticle, ccake::thermo_info::dwds) = particles[iparticle].thermo.dwds;
    host_thermo(iparticle, ccake::thermo_info::dwdB) = particles[iparticle].thermo.dwdB;
    host_thermo(iparticle, ccake::thermo_info::dwdS) = particles[iparticle].thermo.dwdS;
    host_thermo(iparticle, ccake::thermo_info::dwdQ) = particles[iparticle].thermo.dwdQ;
    #ifdef DEBUG_SLOW
    //Couting to check values
    if (iparticle > 998 and iparticle < 1002) {
      std::cout << " t " << particles[iparticle].hydro.t << " gamma " << particles[iparticle].hydro.gamma << \
              " sigma_lab " << particles[iparticle].hydro.sigma_lab << " T " << particles[iparticle].thermo.T \
	     << " thermo.e " << particles[iparticle].thermo.e << " thermo.s " << particles[iparticle].thermo.s << std::endl;
    }
    #endif

  // Copy hydro spacetime matrix
    for (int i = 0; i < 4; ++i){
        for (int j = 0; j < 4; ++j){
            host_hydro_spacetime_matrix(iparticle, ccake::hydro_info::shv, i, j) = particles[iparticle].hydro.shv(i, j);
            host_hydro_spacetime_matrix(iparticle, ccake::hydro_info::sigma_tensor, i, j) = particles[iparticle].hydro.sigma_tensor(i, j);
      }
    }

    // Copy hydro shear auxiliary vectors and matrices
    for (int i = 0; i < 2; ++i) {
        for (int j = i; j < 3; ++j) {
            host_hydro_shear_aux_vector(iparticle, ccake::hydro_info::F_extensive_shear, i, j) = particles[iparticle].hydro.F_extensive_shear(i, j);
            host_hydro_shear_aux_vector(iparticle, ccake::hydro_info::F_sigma_tensor, i, j) = particles[iparticle].hydro.F_sigma_tensor(i, j);
            host_hydro_shear_aux_vector(iparticle, ccake::hydro_info::d_extensive_shv_dt, i, j) = particles[iparticle].hydro.d_extensive_shv_dt(i, j);
            host_hydro_shear_aux_vector(iparticle, ccake::hydro_info::extensive_shv, i, j) = particles[iparticle].hydro.extensive_shv(i, j);
            for (int k = 0; k < D; ++k) {
                int linear_index = j * D + i;
                host_hydro_shear_aux_matrix(iparticle, ccake::hydro_info::M_extensive_shear, i, linear_index) = particles[iparticle].hydro.M_extensive_shear(i, linear_index);
                host_hydro_shear_aux_matrix(iparticle, ccake::hydro_info::M_sigma_tensor, i, linear_index) = particles[iparticle].hydro.M_sigma_tensor(i, linear_index);
            }
            for (int k = 0; k < 3; ++k) {
                int linear_index = j * 3 + i;
                host_hydro_shear_aux_matrix(iparticle, ccake::hydro_info::R_extensive_shear, i, linear_index) = particles[iparticle].hydro.R_extensive_shear(i, linear_index);
            }
        }
    }

    host_input(iparticle, ccake::densities_info::s) = particles[iparticle].input.s;
    host_input(iparticle, ccake::densities_info::rhoB) = particles[iparticle].input.rhoB;
    host_input(iparticle, ccake::densities_info::rhoS) = particles[iparticle].input.rhoS;
    host_input(iparticle, ccake::densities_info::rhoQ) = particles[iparticle].input.rhoQ;

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

    host_extensive(iparticle, ccake::densities_info::s) = particles[iparticle].extensive.s;
    host_extensive(iparticle, ccake::densities_info::rhoB) = particles[iparticle].extensive.rhoB;
    host_extensive(iparticle, ccake::densities_info::rhoS) = particles[iparticle].extensive.rhoS;
    host_extensive(iparticle, ccake::densities_info::rhoQ) = particles[iparticle].extensive.rhoQ;

    host_d_dt_extensive(iparticle, ccake::densities_info::s) = particles[iparticle].d_dt_extensive.s;
    host_d_dt_extensive(iparticle, ccake::densities_info::rhoB) = particles[iparticle].d_dt_extensive.rhoB;
    host_d_dt_extensive(iparticle, ccake::densities_info::rhoS) = particles[iparticle].d_dt_extensive.rhoS;
    host_d_dt_extensive(iparticle, ccake::densities_info::rhoQ) = particles[iparticle].d_dt_extensive.rhoQ;

    host_sph_mass(iparticle, ccake::densities_info::s) = particles[iparticle].sph_mass.s;
    host_sph_mass(iparticle, ccake::densities_info::rhoB) = particles[iparticle].sph_mass.rhoB;
    host_sph_mass(iparticle, ccake::densities_info::rhoS) = particles[iparticle].sph_mass.rhoS;
    host_sph_mass(iparticle, ccake::densities_info::rhoQ) = particles[iparticle].sph_mass.rhoQ;

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
      particles[id].hydro.gradV(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::gradV, i, j);
      particles[id].hydro.M_u(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::M_u, i, j);
      particles[id].hydro.R_u(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::R_u, i, j);
      particles[id].hydro.R_0i_shear(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::R_0i_shear, i, j);
      particles[id].hydro.M_0i_shear(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::M_0i_shear, i, j);
    }
    particles[id].hydro.t = host_hydro_scalar(iparticle, ccake::hydro_info::t);
    particles[id].hydro.causality = host_hydro_scalar(iparticle, ccake::hydro_info::causality);
    particles[id].hydro.bulk = host_hydro_scalar(iparticle, ccake::hydro_info::bulk);
    particles[id].hydro.extensive_bulk = host_hydro_scalar(iparticle, ccake::hydro_info::extensive_bulk);
    particles[id].hydro.a = host_hydro_scalar(iparticle, ccake::hydro_info::a);
    particles[id].hydro.rho_B_ext = host_hydro_scalar(iparticle, ccake::hydro_info::rho_B_ext);
    particles[id].hydro.rho_S_ext = host_hydro_scalar(iparticle, ccake::hydro_info::rho_S_ext);
    particles[id].hydro.rho_Q_ext = host_hydro_scalar(iparticle, ccake::hydro_info::rho_Q_ext);
    particles[id].hydro.tau_Pi = host_hydro_scalar(iparticle, ccake::hydro_info::tau_Pi);
    particles[id].hydro.tau_pi = host_hydro_scalar(iparticle, ccake::hydro_info::tau_pi);
    particles[id].hydro.zeta_Pi = host_hydro_scalar(iparticle, ccake::hydro_info::zeta_Pi);
    particles[id].hydro.sigma_lab = host_hydro_scalar(iparticle, ccake::hydro_info::sigma_lab);
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
    particles[id].hydro.tmunu_trace = host_hydro_scalar(iparticle, ccake::hydro_info::tmunu_trace);
    particles[id].hydro.d_extensive_bulk_dt = host_hydro_scalar(iparticle, ccake::hydro_info::d_extensive_bulk_dt);
    particles[id].hydro.F_extensive_bulk = host_hydro_scalar(iparticle, ccake::hydro_info::F_extensive_bulk);
    particles[id].hydro.F_shv_nabla_u = host_hydro_scalar(iparticle, ccake::hydro_info::F_shv_nabla_u);
    particles[id].hydro.F_extensive_entropy = host_hydro_scalar(iparticle, ccake::hydro_info::F_extensive_entropy);
    particles[id].hydro.shv_nabla_u = host_hydro_scalar(iparticle, ccake::hydro_info::shv_nabla_u);
    particles[id].hydro.j0_ext = host_hydro_scalar(iparticle, ccake::hydro_info::j0_ext);
    particles[id].hydro.shear_knudsen = host_hydro_scalar(iparticle, ccake::hydro_info::shear_knudsen);
    particles[id].hydro.bulk_knudsen = host_hydro_scalar(iparticle, ccake::hydro_info::bulk_knudsen);
    particles[id].hydro.inverse_reynolds_shear = host_hydro_scalar(iparticle, ccake::hydro_info::inverse_reynolds_shear);
    particles[id].hydro.inverse_reynolds_bulk = host_hydro_scalar(iparticle, ccake::hydro_info::inverse_reynolds_bulk);

    for (int i=0; i<D; i++){
      particles[id].hydro.v(i) = host_hydro_vector(iparticle, ccake::hydro_info::v,i);
      particles[id].hydro.u(i) = host_hydro_vector(iparticle, ccake::hydro_info::u,i);
      particles[id].hydro.gradP(i) = host_hydro_vector(iparticle, ccake::hydro_info::gradP,i);
      particles[id].hydro.gradE(i) = host_hydro_vector(iparticle, ccake::hydro_info::gradE,i);
      particles[id].hydro.gradBulk(i) = host_hydro_vector(iparticle, ccake::hydro_info::gradBulk,i);
      particles[id].hydro.divshear(i) = host_hydro_vector(iparticle, ccake::hydro_info::divshear,i);
      particles[id].hydro.gradshear(i) = host_hydro_vector(iparticle, ccake::hydro_info::gradshear,i);
      particles[id].hydro.du_dt(i) = host_hydro_vector(iparticle, ccake::hydro_info::du_dt,i);
      particles[id].hydro.M_extensive_bulk(i) = host_hydro_vector(iparticle, ccake::hydro_info::M_extensive_bulk,i);
      particles[id].hydro.M_shv_nabla_u(i) = host_hydro_vector(iparticle, ccake::hydro_info::M_shv_nabla_u,i);
      particles[id].hydro.M_extensive_entropy(i) = host_hydro_vector(iparticle, ccake::hydro_info::M_extensive_entropy,i);
      particles[id].hydro.F_u(i) = host_hydro_vector(iparticle, ccake::hydro_info::F_u,i);
      particles[id].hydro.F_0i_shear(i) = host_hydro_vector(iparticle, ccake::hydro_info::F_0i_shear,i);
      particles[id].hydro.j_ext(i) = host_hydro_vector(iparticle, ccake::hydro_info::j_ext,i);
      particles[id].r(i) = host_position(iparticle, i);
    }

    //charges
    for(int i=0; i<3; ++i){
      particles[id].hydro.F_extensive_N(i) = host_hydro_vector(iparticle, ccake::hydro_info::F_extensive_N,i);
      particles[id].hydro.R_extensive_entropy(i) = host_hydro_vector(iparticle, ccake::hydro_info::R_extensive_entropy,i);
      particles[id].hydro.R_extensive_bulk(i) = host_hydro_vector(iparticle, ccake::hydro_info::R_extensive_bulk,i);
      for(int j=0; j<3; ++j){
        particles[id].hydro.R_extensive_N(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::R_extensive_N,i,j);
      }
      for(int idir=0; idir<D; idir++){
        particles[id].hydro.M_extensive_N(idir,i) = host_hydro_space_matrix(iparticle, ccake::hydro_info::M_extensive_N,idir,i);
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
    particles[id].thermo.dwds = host_thermo(iparticle, ccake::thermo_info::dwds);
    particles[id].thermo.dwdB = host_thermo(iparticle, ccake::thermo_info::dwdB);
    particles[id].thermo.dwdS = host_thermo(iparticle, ccake::thermo_info::dwdS);
    particles[id].thermo.dwdQ = host_thermo(iparticle, ccake::thermo_info::dwdQ);

    std::cout << "Made it to " << __PRETTY_FUNCTION__ << ":: line = " << __LINE__ << std::endl;
    for (int i=0; i<4; i++){
    for (int j=0; j<4; j++){
      particles[id].hydro.shv(i,j) = host_hydro_spacetime_matrix(iparticle, ccake::hydro_info::shv, i, j);
      particles[id].hydro.sigma_tensor(i,j) = host_hydro_spacetime_matrix(iparticle, ccake::hydro_info::sigma_tensor, i, j);
    }}
    for (int i=0; i<2; i++){
      for (int j=i; j<3; j++){
        particles[id].hydro.F_extensive_shear(i,j) = host_hydro_shear_aux_vector(iparticle, ccake::hydro_info::F_extensive_shear, i, j);
        particles[id].hydro.F_sigma_tensor(i,j) = host_hydro_shear_aux_vector(iparticle, ccake::hydro_info::F_sigma_tensor, i, j);
        particles[id].hydro.d_extensive_shv_dt(i,j) = host_hydro_shear_aux_vector(iparticle, ccake::hydro_info::d_extensive_shv_dt, i, j);
        particles[id].hydro.extensive_shv(i,j) = host_hydro_shear_aux_vector(iparticle, ccake::hydro_info::extensive_shv, i, j);
        for (int k=0; k<D; k++){
          int linear_index = j * D + k;
          particles[id].hydro.M_sigma_tensor(i,linear_index) = host_hydro_shear_aux_matrix(iparticle, ccake::hydro_info::M_sigma_tensor, i, linear_index);
          particles[id].hydro.M_extensive_shear(i,linear_index) = host_hydro_shear_aux_matrix(iparticle, ccake::hydro_info::M_extensive_shear, i, linear_index);
        }
        for (int k=0; k<3; k++){
          int linear_index = j * 3 + k;
          particles[id].hydro.R_extensive_shear(i,linear_index) = host_hydro_shear_aux_matrix(iparticle, ccake::hydro_info::R_extensive_shear, i, linear_index);
        }
      }
    }


    std::cout << "Made it to " << __PRETTY_FUNCTION__ << ":: line = " << __LINE__ << std::endl;

    particles[id].input.s = host_input(iparticle, ccake::densities_info::s);
    particles[id].input.rhoB = host_input(iparticle, ccake::densities_info::rhoB);
    particles[id].input.rhoS = host_input(iparticle, ccake::densities_info::rhoS);
    particles[id].input.rhoQ = host_input(iparticle, ccake::densities_info::rhoQ);

    std::cout << "Made it to " << __PRETTY_FUNCTION__ << ":: line = " << __LINE__ << std::endl;
    particles[id].smoothed.s = host_smoothed(iparticle, ccake::densities_info::s);
    particles[id].smoothed.rhoB = host_smoothed(iparticle, ccake::densities_info::rhoB);
    particles[id].smoothed.rhoS = host_smoothed(iparticle, ccake::densities_info::rhoS);
    particles[id].smoothed.rhoQ = host_smoothed(iparticle, ccake::densities_info::rhoQ);

    std::cout << "Made it to " << __PRETTY_FUNCTION__ << ":: line = " << __LINE__ << std::endl;
    particles[id].extensive.s = host_extensive(iparticle, ccake::densities_info::s);
    particles[id].extensive.rhoB = host_extensive(iparticle, ccake::densities_info::rhoB);
    particles[id].extensive.rhoS = host_extensive(iparticle, ccake::densities_info::rhoS);
    particles[id].extensive.rhoQ = host_extensive(iparticle, ccake::densities_info::rhoQ);

    std::cout << "Made it to " << __PRETTY_FUNCTION__ << ":: line = " << __LINE__ << std::endl;
    particles[id].d_dt_extensive.s = host_d_dt_extensive(iparticle, ccake::densities_info::s);
    particles[id].d_dt_extensive.rhoB = host_d_dt_extensive(iparticle, ccake::densities_info::rhoB);
    particles[id].d_dt_extensive.rhoS = host_d_dt_extensive(iparticle, ccake::densities_info::rhoS);
    particles[id].d_dt_extensive.rhoQ = host_d_dt_extensive(iparticle, ccake::densities_info::rhoQ);

    std::cout << "Made it to " << __PRETTY_FUNCTION__ << ":: line = " << __LINE__ << std::endl;
    particles[id].sph_mass.s = host_sph_mass(iparticle, ccake::densities_info::s);
    particles[id].sph_mass.rhoB = host_sph_mass(iparticle, ccake::densities_info::rhoB);
    particles[id].sph_mass.rhoS = host_sph_mass(iparticle, ccake::densities_info::rhoS);
    particles[id].sph_mass.rhoQ = host_sph_mass(iparticle, ccake::densities_info::rhoQ);

    std::cout << "Made it to " << __PRETTY_FUNCTION__ << ":: line = " << __LINE__ << std::endl;
    particles[id].efcheck = host_efcheck(iparticle);
    particles[id].contribution_to_total_E = host_contribution_to_total_E(iparticle);
    particles[id].contribution_to_total_dEz = host_contribution_to_total_dEz(iparticle);
    particles[id].contribution_to_total_Ez = host_contribution_to_total_Ez(iparticle);
    particles[id].ID = host_id(iparticle);
    particles[id].btrack = host_btrack(iparticle);
    particles[id].Freeze = host_freeze(iparticle);
  }
  std::cout << "Made it to " << __PRETTY_FUNCTION__ << ":: line = " << __LINE__ << std::endl;

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
    min_pos[idir] *= 2.; //Grid must be 100% extensiveger ///TODO: Allow this to be an optional input parameter

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
  /// Maybe this paralle for fix it?
  ///Need a way to call Cabana::NeighborList::numNeighbor directly from GPU
  ///without seg fault
  Kokkos::parallel_for("UpdateNeighbors", Kokkos::RangePolicy<>(0, n_particles), KOKKOS_LAMBDA(int i) {
    device_btrack(i) = (device_btrack(i) == -1) ? -1
                      : Cabana::NeighborList<ListType>::numNeighbor(neighbour_list, i);
  });
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
    local_S += device_extensive(i, ccake::densities_info::s)*device_sph_mass(i, ccake::densities_info::s);
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
    Btotal += device_extensive(i, ccake::densities_info::rhoB)*device_sph_mass(i, ccake::densities_info::rhoB);
  };
  auto get_total_S = KOKKOS_LAMBDA(const int &i, double &Stotal)
  {
    Stotal += device_extensive(i, ccake::densities_info::rhoS)*device_sph_mass(i, ccake::densities_info::rhoS);
  };
  auto get_total_Q = KOKKOS_LAMBDA(const int &i, double &Qtotal)
  {
    Qtotal += device_extensive(i, ccake::densities_info::rhoQ)*device_sph_mass(i, ccake::densities_info::rhoQ);
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
void SystemState<D>::conservation_energy(bool first_iteration, double t)
{
  ///////////////////////////////////////////////
  // don't bother checking energy conservation on
  // intermediate RK steps
  // calculate total energy (T^{00})

  CREATE_VIEW(device_, cabana_particles);
  E = 0.0;
  auto get_total_energy = KOKKOS_LAMBDA(const int i, double &local_E)
  {
    double w = device_thermo(i, ccake::thermo_info::w);
    double dE = device_contribution_to_total_E(i);
    double p = device_thermo(i, ccake::thermo_info::p);
    double sph_mass = device_sph_mass(i, ccake::densities_info::s);
    double sigma_lab = device_hydro_scalar(i, ccake::hydro_info::sigma_lab);
    double bulk = device_hydro_scalar(i, ccake::hydro_info::bulk);
    double shv00 = device_hydro_spacetime_matrix(i, ccake::hydro_info::shv, 0,0);
    double gamma = device_hydro_scalar(i, ccake::hydro_info::gamma);
    double sigma = device_hydro_scalar(i, ccake::hydro_info::sigma);
    double g2 = gamma*gamma;

    double C = w + bulk;

    local_E += (C * g2 - p - bulk + shv00) * sph_mass * t / sigma_lab;
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
