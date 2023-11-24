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


////////////////////////////////////////////////////////////////////////////////
template<unsigned int D>
void SystemState<D>::initialize()  // formerly called "manualenter"
{
  formatted_output::report("Initializing system");

  t = settingsPtr->t0;
  hT = settingsPtr->hT;

  formatted_output::update("set freeze out parameters");

  formatted_output::detail("freeze out temperature = "
                           + to_string(settingsPtr->Freeze_Out_Temperature
                                        *hbarc_MeVfm) + " MeV");
  return;
}

template<unsigned int D>
void SystemState<D>::allocate_cabana_particles(){
    formatted_output::detail("Initializing device memory");
    ///Allocate memory for the particles
    cabana_particles = Cabana::AoSoA<CabanaParticle, DeviceType, VECTOR_LENGTH>("particles", n_particles);
    copy_host_to_device();
}

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
      host_hydro_space_matrix(iparticle, ccake::hydro_info::Imat, i, j) = particles[iparticle].hydro.Imat(i,j);
      host_hydro_space_matrix(iparticle, ccake::hydro_info::gradV, i, j) = particles[iparticle].hydro.gradV(i,j);
      host_hydro_space_matrix(iparticle, ccake::hydro_info::gradU, i, j) = particles[iparticle].hydro.gradU(i,j);
      host_hydro_space_matrix(iparticle, ccake::hydro_info::uu, i, j) = particles[iparticle].hydro.uu(i,j);
      host_hydro_space_matrix(iparticle, ccake::hydro_info::pimin, i, j) = particles[iparticle].hydro.pimin(i,j);
      host_hydro_space_matrix(iparticle, ccake::hydro_info::piu, i, j) = particles[iparticle].hydro.piu(i,j);
      host_hydro_space_matrix(iparticle, ccake::hydro_info::pimin, i, j) = particles[iparticle].hydro.pimin(i,j);
      host_hydro_space_matrix(iparticle, ccake::hydro_info::piutot, i, j) = particles[iparticle].hydro.piutot(i,j);
      host_hydro_space_matrix(iparticle, ccake::hydro_info::shv1, i, j) = particles[iparticle].hydro.shv1(i,j);
      host_hydro_space_matrix(iparticle, ccake::hydro_info::shv2, i, j) = particles[iparticle].hydro.shv2(i,j);
      host_hydro_space_matrix(iparticle, ccake::hydro_info::shv3, i, j) = particles[iparticle].hydro.shv3(i,j);
      host_hydro_space_matrix(iparticle, ccake::hydro_info::shv4, i, j) = particles[iparticle].hydro.shv4(i,j);
      host_hydro_space_matrix(iparticle, ccake::hydro_info::dshv_dt, i, j) = particles[iparticle].hydro.dshv_dt(i,j);
    }
    host_hydro_scalar(iparticle, ccake::hydro_info::t) = particles[iparticle].hydro.t;
    host_hydro_scalar(iparticle, ccake::hydro_info::Agam) = particles[iparticle].hydro.Agam;
    host_hydro_scalar(iparticle, ccake::hydro_info::Agam2) = particles[iparticle].hydro.Agam2;
    host_hydro_scalar(iparticle, ccake::hydro_info::shv33) = particles[iparticle].hydro.shv33;
    host_hydro_scalar(iparticle, ccake::hydro_info::gamma) = particles[iparticle].hydro.gamma;
    host_hydro_scalar(iparticle, ccake::hydro_info::Bulk) = particles[iparticle].hydro.Bulk;
    host_hydro_scalar(iparticle, ccake::hydro_info::bigPI) = particles[iparticle].hydro.bigPI;
    host_hydro_scalar(iparticle, ccake::hydro_info::C) = particles[iparticle].hydro.C;
    host_hydro_scalar(iparticle, ccake::hydro_info::tauRelax) = particles[iparticle].hydro.tauRelax;
    host_hydro_scalar(iparticle, ccake::hydro_info::stauRelax) = particles[iparticle].hydro.stauRelax;
    host_hydro_scalar(iparticle, ccake::hydro_info::zeta) = particles[iparticle].hydro.zeta;
    host_hydro_scalar(iparticle, ccake::hydro_info::setas) = particles[iparticle].hydro.setas;
    host_hydro_scalar(iparticle, ccake::hydro_info::Ctot) = particles[iparticle].hydro.Ctot;
    host_hydro_scalar(iparticle, ccake::hydro_info::Btot) = particles[iparticle].hydro.Btot;
    host_hydro_scalar(iparticle, ccake::hydro_info::sigma) = particles[iparticle].hydro.sigma;
    host_hydro_scalar(iparticle, ccake::hydro_info::dsigma_dt) = particles[iparticle].hydro.dsigma_dt;
    host_hydro_scalar(iparticle, ccake::hydro_info::gamma_squared) = particles[iparticle].hydro.gamma_squared;
    host_hydro_scalar(iparticle, ccake::hydro_info::gamma_cube) = particles[iparticle].hydro.gamma_cube;
    host_hydro_scalar(iparticle, ccake::hydro_info::gamma_tau) = particles[iparticle].hydro.gamma_tau;
    host_hydro_scalar(iparticle, ccake::hydro_info::eta_o_tau) = particles[iparticle].hydro.eta_o_tau;
    host_hydro_scalar(iparticle, ccake::hydro_info::dwdsT1) = particles[iparticle].hydro.dwdsT1;
    host_hydro_scalar(iparticle, ccake::hydro_info::sigl) = particles[iparticle].hydro.sigl;
    host_hydro_scalar(iparticle, ccake::hydro_info::varsigma) = particles[iparticle].hydro.varsigma;
    host_hydro_scalar(iparticle, ccake::hydro_info::bigtheta) = particles[iparticle].hydro.bigtheta;
    host_hydro_scalar(iparticle, ccake::hydro_info::inside) = particles[iparticle].hydro.inside;
    host_hydro_scalar(iparticle, ccake::hydro_info::div_u) = particles[iparticle].hydro.div_u;
    host_hydro_scalar(iparticle, ccake::hydro_info::dBulk_dt) = particles[iparticle].hydro.dBulk_dt;
    for (int i=0; i<D; i++){
      host_hydro_vector(iparticle, ccake::hydro_info::v,i) = particles[iparticle].hydro.v(i);
      host_hydro_vector(iparticle, ccake::hydro_info::u,i) = particles[iparticle].hydro.u(i);
      host_hydro_vector(iparticle, ccake::hydro_info::gradP,i) = particles[iparticle].hydro.gradP(i);
      host_hydro_vector(iparticle, ccake::hydro_info::gradE,i) = particles[iparticle].hydro.gradE(i);
      host_hydro_vector(iparticle, ccake::hydro_info::gradBulk,i) = particles[iparticle].hydro.gradBulk(i);
      host_hydro_vector(iparticle, ccake::hydro_info::divshear,i) = particles[iparticle].hydro.divshear(i);
      host_hydro_vector(iparticle, ccake::hydro_info::gradshear,i) = particles[iparticle].hydro.gradshear(i);
      host_hydro_vector(iparticle, ccake::hydro_info::du_dt,i) = particles[iparticle].hydro.du_dt(i);

      host_position(iparticle, i) = particles[iparticle].r(i);
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

    for (int i=0; i<D+1; i++)
    for (int j=0; j<D+1; j++)
      host_hydro_spacetime_matrix(iparticle, ccake::hydro_info::shv, i, j) =  particles[iparticle].hydro.shv(i,j);

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

/// @brief Takes the data out of the device memory (managed by cabana) and
/// back to particle data structure
/// @details This function is expensive computationally  and its use should
/// be avoided if possible. It creates a temporary AoSoA in host memory space
/// and mirrors the device memory. Just then it populates the particles structure
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
      particles[id].hydro.gradU(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::gradU, i, j);
      particles[id].hydro.uu(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::uu, i, j);
      particles[id].hydro.pimin(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::pimin, i, j);
      particles[id].hydro.piu(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::piu, i, j);
      particles[id].hydro.pimin(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::pimin, i, j);
      particles[id].hydro.piutot(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::piutot, i, j);
      particles[id].hydro.shv1(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::shv1, i, j);
      particles[id].hydro.shv2(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::shv2, i, j);
      particles[id].hydro.shv3(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::shv3, i, j);
      particles[id].hydro.shv4(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::shv4, i, j);
      particles[id].hydro.dshv_dt(i,j) = host_hydro_space_matrix(iparticle, ccake::hydro_info::dshv_dt, i, j);
    }
    particles[id].hydro.t = host_hydro_scalar(iparticle, ccake::hydro_info::t);
    particles[id].hydro.Agam = host_hydro_scalar(iparticle, ccake::hydro_info::Agam);
    particles[id].hydro.Agam2 = host_hydro_scalar(iparticle, ccake::hydro_info::Agam2);
    particles[id].hydro.shv33 = host_hydro_scalar(iparticle, ccake::hydro_info::shv33);
    particles[id].hydro.gamma = host_hydro_scalar(iparticle, ccake::hydro_info::gamma);
    particles[id].hydro.Bulk = host_hydro_scalar(iparticle, ccake::hydro_info::Bulk);
    particles[id].hydro.bigPI = host_hydro_scalar(iparticle, ccake::hydro_info::bigPI);
    particles[id].hydro.C = host_hydro_scalar(iparticle, ccake::hydro_info::C);
    particles[id].hydro.tauRelax = host_hydro_scalar(iparticle, ccake::hydro_info::tauRelax);
    particles[id].hydro.stauRelax = host_hydro_scalar(iparticle, ccake::hydro_info::stauRelax);
    particles[id].hydro.zeta = host_hydro_scalar(iparticle, ccake::hydro_info::zeta);
    particles[id].hydro.setas = host_hydro_scalar(iparticle, ccake::hydro_info::setas);
    particles[id].hydro.Ctot = host_hydro_scalar(iparticle, ccake::hydro_info::Ctot);
    particles[id].hydro.Btot = host_hydro_scalar(iparticle, ccake::hydro_info::Btot);
    particles[id].hydro.sigma = host_hydro_scalar(iparticle, ccake::hydro_info::sigma);
    particles[id].hydro.dsigma_dt = host_hydro_scalar(iparticle, ccake::hydro_info::dsigma_dt);
    particles[id].hydro.gamma_squared = host_hydro_scalar(iparticle, ccake::hydro_info::gamma_squared);
    particles[id].hydro.gamma_cube = host_hydro_scalar(iparticle, ccake::hydro_info::gamma_cube);
    particles[id].hydro.gamma_tau = host_hydro_scalar(iparticle, ccake::hydro_info::gamma_tau);
    particles[id].hydro.eta_o_tau = host_hydro_scalar(iparticle, ccake::hydro_info::eta_o_tau);
    particles[id].hydro.dwdsT1 = host_hydro_scalar(iparticle, ccake::hydro_info::dwdsT1);
    particles[id].hydro.sigl = host_hydro_scalar(iparticle, ccake::hydro_info::sigl);
    particles[id].hydro.varsigma = host_hydro_scalar(iparticle, ccake::hydro_info::varsigma);
    particles[id].hydro.bigtheta = host_hydro_scalar(iparticle, ccake::hydro_info::bigtheta);
    particles[id].hydro.inside = host_hydro_scalar(iparticle, ccake::hydro_info::inside);
    particles[id].hydro.div_u = host_hydro_scalar(iparticle, ccake::hydro_info::div_u);
    particles[id].hydro.dBulk_dt = host_hydro_scalar(iparticle, ccake::hydro_info::dBulk_dt);
    for (int i=0; i<D; i++){
      particles[id].hydro.v(i) = host_hydro_vector(iparticle, ccake::hydro_info::v,i);
      particles[id].hydro.u(i) = host_hydro_vector(iparticle, ccake::hydro_info::u,i);
      particles[id].hydro.gradP(i) = host_hydro_vector(iparticle, ccake::hydro_info::gradP,i);
      particles[id].hydro.gradE(i) = host_hydro_vector(iparticle, ccake::hydro_info::gradE,i);
      particles[id].hydro.gradBulk(i) = host_hydro_vector(iparticle, ccake::hydro_info::gradBulk,i);
      particles[id].hydro.divshear(i) = host_hydro_vector(iparticle, ccake::hydro_info::divshear,i);
      particles[id].hydro.gradshear(i) = host_hydro_vector(iparticle, ccake::hydro_info::gradshear,i);
      particles[id].hydro.du_dt(i) = host_hydro_vector(iparticle, ccake::hydro_info::du_dt,i);

      particles[id].r(i) = host_position(iparticle, i);
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

    for (int i=0; i<D+1; i++)
    for (int j=0; j<D+1; j++)
      particles[id].hydro.shv(i,j) = host_hydro_spacetime_matrix(iparticle, ccake::hydro_info::shv, i, j);

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
/// @details This function should be called after the particles have been moved.
/// Beware that this operation is very expensive and calls to this function should
/// be minimized. It is also a good candidate for optimization. In a 11th Gen
/// Intel Core i7-11370H running under WSL 2, using the apptainer images, the
/// inclusion of single call to this function in the main loop increased the
/// execution time by ~1s (in a trento event generated with
/// `trento --grid-max 10 --grid-step .1 Pb Pb --b-max=3 --random-seed=42`)
/// @tparam D
template<unsigned int D>
void SystemState<D>::reset_neighbour_list(){
  double min_pos[3], max_pos[3], delta[3];
  min_pos[0] = settingsPtr->xmin;
  min_pos[1] = settingsPtr->ymin;
  min_pos[2] = settingsPtr->etamin;
  for(int idir=0; idir<3; ++idir)
    min_pos[idir] *= 2.; //Grid must be 100% bigger ///TODO: Allow this to be an optional input parameter
  for(int idir=0; idir<3; ++idir)
    max_pos[idir] = -min_pos[idir];
  for(int idir=0; idir<3; ++idir)
    delta[idir] = 2.*settingsPtr->hT;

  CREATE_VIEW(device_, cabana_particles);

  //Cabana::LinkedCellList<DeviceType> cell_list(device_position, delta, min_pos, max_pos);
  //Cabana::permute( cell_list, cabana_particles ); ///TODO: Enabling may result in potential performance gain. Needs testing.
  double neighborhood_radius = 2*settingsPtr->hT;
  double cell_ratio = 1.; //neighbour to cell_space ratio
  neighbour_list = ListType (  device_position, 0, device_position.size(),
                                          neighborhood_radius, cell_ratio, min_pos, max_pos
                           );
  Kokkos::fence();
}

template<unsigned int D>
void SystemState<D>::initialize_linklist()
{
  formatted_output::report("Initializing linklist");

  reset_neighbour_list();
  return;
}

///////////////////////////////////////
//TODO: Parallelize with Kokkos::parallel_reduce
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

///////////////////////////////////////
//TODO: Parallelize with Kokkos::parallel_reduce
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
  if ( rk2 == 1 )
  {
    // calculate total energy (T^{00})
    E = 0.0;
    for ( auto & p : particles )
    {
      p.contribution_to_total_E
         = ( p.hydro.C*p.hydro.gamma_squared - p.p() - p.hydro.bigPI + p.hydro.shv(0,0) )
           * p.norm_spec.s * t / p.hydro.sigma;
      E += p.contribution_to_total_E;
    }

    // store initial total energy
    // for checking subsequent energy loss
    //if (linklist.first==1)
    {
      //linklist.first = 0;
      E0             = E;
    }

    // Ez is initially set to zero,
    // updated subsequently during RK integration
    Etot  = E + Ez;
    Eloss = (E0-Etot)/E0*100;
    rk2   = 0;
  }

  ///////////////////////////////////////////////
  // this enters the RK routine and should be
  // done for intermediate steps as well;
  // this gives the longitudinal energy flux (~T^{\eta\eta})
  dEz = 0.0;
  double t2 = t*t;
  for ( auto & p : particles )
  {
    p.contribution_to_total_dEz
         = ( p.p() + p.hydro.bigPI + p.hydro.shv33*t2 )
           * p.norm_spec.s / p.hydro.sigma;
    dEz += p.contribution_to_total_dEz;
  }

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
    double p_plus_Pi = p.p() + p.hydro.bigPI;
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