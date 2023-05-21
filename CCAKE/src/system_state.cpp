#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>

#include "system_state.h"

//using namespace std;
using std::cout;
using std::endl;
using std::string;

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

  ///Allocate memory for the particles
  cabana_particles = Cabana::AoSoA<CabanaParticle, DeviceType, VECTOR_LENGTH>("particles", n_particles);
  using SerialHost = Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>;
  //Auxiliary AoSoA for copying particles from/to host
  Cabana::AoSoA<CabanaParticle, SerialHost, 8> particles_h("particles_h",n_particles);
  formatted_output::detail("Initializing device memory");
  //Create vies to particles
  CREATE_VIEW(host_, particles_h);
  CREATE_VIEW(device_, cabana_particles);

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
    host_hydro_scalar(iparticle, ccake::hydro_info::g2) = particles[iparticle].hydro.g2;
    host_hydro_scalar(iparticle, ccake::hydro_info::g3) = particles[iparticle].hydro.g3;
    host_hydro_scalar(iparticle, ccake::hydro_info::gt) = particles[iparticle].hydro.gt;
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
      host_hydro_vector(iparticle, ccake::hydro_info::gradBulk,i) = particles[iparticle].hydro.gradBulk(i);
      host_hydro_vector(iparticle, ccake::hydro_info::divshear,i) = particles[iparticle].hydro.divshear(i);
      host_hydro_vector(iparticle, ccake::hydro_info::gradshear,i) = particles[iparticle].hydro.gradshear(i);
      host_hydro_vector(iparticle, ccake::hydro_info::du_dt,i) = particles[iparticle].hydro.du_dt(i);

      host_position(iparticle, i) = particles[iparticle].r(i);
    }
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

}

template<unsigned int D>
void SystemState<D>::initialize_linklist()
{
  formatted_output::report("Initializing linklist");
  double min_pos[D], max_pos[D];
  CREATE_VIEW(device_, cabana_particles);

  cout << "Number of particles: " << n_particles << endl;
  cout << "min x:.............: " << settingsPtr->xmin << endl;
  cout << "min y:.............: " << settingsPtr->ymin << endl;

  //Find max position and minimum position
  double min_aux = 0;
  double max_aux = 0;
  for (int idir=0; idir<D; ++idir){
    Kokkos::parallel_reduce(n_particles,
    KOKKOS_CLASS_LAMBDA (const int& iparticle, double& min_aux, double& max_aux ) {
      double x = device_position(iparticle,0);
      double y = device_position(iparticle,1);
      double eps = device_thermo(iparticle, ccake::thermo_info::e);
      double p = device_position(iparticle,idir);
      //cout << p << endl;
      min_aux = min_aux >= p ? p : min_aux;
      max_aux = max_aux <  p ? p : max_aux;
    }
    , min_aux, max_aux);
    Kokkos::fence();
    min_pos[idir] = min_aux;
    max_pos[idir] = max_aux;
  }

  for(int idir=0; idir<D; ++idir){
    min_pos[idir] *= 2;
    max_pos[idir] *= 2;
  }
  for(int idir=D; idir<3; ++idir){
    min_pos[idir] *= 0;
    max_pos[idir] *= 0;
  }

  double grid_spacing[3] = {2*settingsPtr->hT, 2*settingsPtr->hT, 2*settingsPtr->hEta};
  grid = Cabana::LinkedCellList<DeviceType>(device_position, grid_spacing,
                                             min_pos, max_pos);
  Cabana::permute( grid, cabana_particles );
  Kokkos::fence();

  cout << "Grid edges are { ";
  for(int idir=0; idir<D; ++idir){
    cout << "("<< min_pos[idir] <<", "<<max_pos[idir] << ") ";
  }
  cout << "}" << endl;

  std::ofstream fs("validate_ic.txt");
  for (int ipart=0; ipart<n_particles; ++ipart){
    fs << device_position(ipart,0) << " " << device_position(ipart,1)
         << " " << device_thermo(ipart, ccake::thermo_info::e) << endl;
  }
  fs.close();
  return;
}

///////////////////////////////////////
//TODO: Parallelize with Kokkos::parallel_reduce
template<unsigned int D>
void SystemState<D>::conservation_entropy()
{
  S = 0.0;
  for ( auto & p : particles )
    S += p.specific.s*p.norm_spec.s;

  //if (linklist.first==1)
  S0 = S;
}

///////////////////////////////////////
//TODO: Parallelize with Kokkos::parallel_reduce
template<unsigned int D>
void SystemState<D>::conservation_BSQ()
{
  // reset
  Btotal = 0.0;
  Stotal = 0.0;
  Qtotal = 0.0;

  // sum
  for ( auto & p : particles )
  {
    Btotal += p.specific.rhoB*p.norm_spec.rhoB;
    Stotal += p.specific.rhoS*p.norm_spec.rhoS;
    Qtotal += p.specific.rhoQ*p.norm_spec.rhoQ;
  }

  // save initial totals
  //if (linklist.first==1)
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
void SystemState<D>::conservation_energy()
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
         = ( p.hydro.C*p.hydro.g2 - p.p() - p.hydro.bigPI + p.hydro.shv(0,0) )
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