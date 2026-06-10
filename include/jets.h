#ifndef BBMG_H
#define BBMG_H

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <omp.h>
#include <random>

#include "constants.h"
#include "kernel.h"
#include "particle.h"
#include "settings.h"
#include "system_state.h"
#include "thermodynamic_info.h"
#include "milne.hpp"

#include <Cabana_Core.hpp>



//using namespace std;

namespace ccake
{
  using jets = Cabana::MemberTypes<double[3], // position
                                  double[3], // velocity
                                  double, // initial x
                                  double, // initial y
                                  double, // initial eta
                                  double, // rho0
                                  double,    // rho
                                  double,    // Temp
                                  double,    // Initial Temp
                                  double,    // phi
                                  double,    // rapidity (jet rapidity)
                                  double,    // line integral
                                  double,    // e_loss (energy lost this step)
                                  double,    // gamma
                                  double,    // vmag
                                  double,    // vang
                                  double,    // flow
                                  int,       // PID
                                  int>;      // Frozen out

    namespace jets_enum{
    enum jets_members
    {
      position,
      velocity,
      x,
      y,
      eta,
      rho0,
      rho,
      T,
      T0,
      phi,
      rapidity,
      energy_lost,
      e_loss,
      gam,
      vmag,
      vang,
      flow,
      PID,
      Frozen,
      //NUM_JET_INFO
    };
  //#define JETS_MEMBERS double[ccake::jets_enum::NUM_JET_INFO]
  }

  using jet_array = Cabana::AoSoA<jets, DeviceType,VECTOR_LENGTH>;
  
  #define jets_VIEW(prefix, jets_aosoa) \
  auto CONCAT(prefix, position) = Cabana::slice<jets_enum::position>(jets_aosoa); \
  auto CONCAT(prefix, velocity) = Cabana::slice<jets_enum::velocity>(jets_aosoa); \
  auto CONCAT(prefix, x) = Cabana::slice<jets_enum::x>(jets_aosoa); \
  auto CONCAT(prefix, y) = Cabana::slice<jets_enum::y>(jets_aosoa); \
  auto CONCAT(prefix, eta) = Cabana::slice<jets_enum::eta>(jets_aosoa); \
  auto CONCAT(prefix, rho0) = Cabana::slice<jets_enum::rho0>(jets_aosoa); \
  auto CONCAT(prefix, rho) = Cabana::slice<jets_enum::rho>(jets_aosoa); \
  auto CONCAT(prefix, T) = Cabana::slice<jets_enum::T>(jets_aosoa); \
  auto CONCAT(prefix, T0) = Cabana::slice<jets_enum::T0>(jets_aosoa); \
  auto CONCAT(prefix, phi) = Cabana::slice<jets_enum::phi>(jets_aosoa); \
  auto CONCAT(prefix, rapidity) = Cabana::slice<jets_enum::rapidity>(jets_aosoa); \
  auto CONCAT(prefix, energy_lost) = Cabana::slice<jets_enum::energy_lost>(jets_aosoa); \
  auto CONCAT(prefix, e_loss) = Cabana::slice<jets_enum::e_loss>(jets_aosoa); \
  auto CONCAT(prefix, gam) = Cabana::slice<jets_enum::gam>(jets_aosoa); \
  auto CONCAT(prefix, vmag) = Cabana::slice<jets_enum::vmag>(jets_aosoa); \
  auto CONCAT(prefix, vang) = Cabana::slice<jets_enum::vang>(jets_aosoa); \
  auto CONCAT(prefix, flow) = Cabana::slice<jets_enum::flow>(jets_aosoa); \
  auto CONCAT(prefix, PID) = Cabana::slice<jets_enum::PID>(jets_aosoa); \
  auto CONCAT(prefix, Frozen) = Cabana::slice<jets_enum::Frozen>(jets_aosoa);


template <unsigned int D>
class BBMG {

private:
    std::shared_ptr<Settings> settingsPtr;
    std::shared_ptr<SystemState<D>> systemPtr;

    int z, a, c, q;
    int num_jets;
    double Freezeout_Temp;
    double area;
    double vjet; // Taken to be c for jets
    double Cg, Cq;
    double gridx, gridy;
    
    double Pfg, Pfq; ///\@todo: Move or eliminate based on need for energy loss equation with energy dependence


    //void inter(field &f); // interpolation
    //KOKKOS_FUNCTION
    //double efluc();
    KOKKOS_FUNCTION
    double get_kappa(double T);


public:

    BBMG(){}
    BBMG( std::shared_ptr<Settings> settingsPtr_in, std::shared_ptr<SystemState<D>> systemPtr_in )
            : settingsPtr(settingsPtr_in),
              systemPtr(systemPtr_in){};
    void copy_host_to_device_BBMG();
    void copy_device_to_host_BBMG();
    void propagate();

    //KOKKOS_FUNCTION
    //double flow(field &f);

    void initial();

    Cabana::AoSoA<jets, DeviceType,VECTOR_LENGTH> jetInfo;

    struct field  // All fetched info for jet initialization and propagation
    {
        int sph, on;
        double rho, rho0, T, v[3];
        double r[3], phi, rapidity, energy_lost, e_loss, x, y, eta;
        int PID, Frozen;
        double gam, vmag, vang, flow;
        double T0;
    };

    vector<field> jetInfo_host; ///\@todo: Change these vectors to the AoSoA's
    vector<field> jetFreezeOut;
};
} //ccake namespace end

using namespace ccake;

template class BBMG<2>;
template class BBMG<3>;

template <unsigned int D>
double BBMG<D>::get_kappa(double T) //This is taken from one plot we found of Barbara's, no reason to believe this should be set in stone. Will test multiple Kappa's
{
  return 2.5*std::exp(-9.22*T) + 0.04;
}

template <unsigned int D>
void BBMG<D>::initial()
{
      srand( time( NULL ) );
  Freezeout_Temp    = 150;// Decoherence temperature in MeV
  Cg                = 3; // Cassimir const gluons
  Cq                = 4./3; // Cassimir const quarks
  // All next quantities are part of the BBMG parameters
  q                 = settingsPtr->jets_Fluctuations; // Fluctuation parameter, subject to change
  z                 = settingsPtr->jets_Length_scaling; // path length dependence
  a                 = settingsPtr->jets_Energy_scaling; //Initial jet energy dependence
  c                 = (2+z-a)/3; //medium temperature dependence
  
  //===============================================
  vjet  = 1;

  gridx = settingsPtr->stepx;
  gridy = settingsPtr->stepy;
  Pfg = 10 / constants::hbarc_GeVfm;
  Pfq = 10 / constants::hbarc_GeVfm;

  auto& p = systemPtr->particles;

  // Push one jet + its back-to-back partner. Caller supplies local T and ρ.
  auto push_pair = [&](double xj, double yj, double etj,
                       double phij, double rapidity_,
                       double rho0, double T_MeV, int pid_base) {
    double line0 = 0.5 * get_kappa(T_MeV / 1000)
                       * std::pow(settingsPtr->t0, z) * std::pow(rho0, c) * settingsPtr->dt;
    for (int j = 0; j < 2; ++j) {
      field f{};
      f.r[0] = xj; f.r[1] = yj;
      if (D == 3) f.r[2] = etj;
      f.x = xj; f.y = yj; f.eta = (D == 3) ? etj : 0.0;
      f.phi   = phij + j * PI;
      f.rapidity = (D == 3) ? rapidity_ : 0.0;
      f.rho0 = rho0; f.T = T_MeV; f.T0 = T_MeV;
      f.energy_lost = line0;
      f.e_loss = 0.0;
      f.PID = 2 * pid_base + j;
      f.Frozen = 1;
      jetInfo_host.push_back(f);
    }
  };

  // Kernel-weighted SPH sum at (xj, yj, etj) — used by file and hardcoded
  // modes where the jet position doesn't coincide with an SPH particle.
  // Returns (rho0, T_MeV); falls back to nearest particle if no neighbors
  // are within the kernel support.
  auto sph_average = [&](double xj, double yj, double etj) {
    double r_jet[D];
    r_jet[0] = xj; r_jet[1] = yj;
    if (D == 3) r_jet[2] = etj;
    double norm = 0.0, T_sum = 0.0, rho_sum = 0.0;
    for (size_t k = 0; k < p.size(); ++k) {
      double r_sph[D];
      for (int d = 0; d < (int)D; ++d) r_sph[d] = p[k].r(d);
      double w = SPHkernel<D>::kernel(r_sph, r_jet,
                                      settingsPtr->hT, settingsPtr->hEta);
      if (w <= 0.0) continue;
      norm    += w;
      T_sum   += p[k].T() * constants::hbarc_MeVfm * w;
      rho_sum += (p[k].p() / p[k].T()) * w;
    }
    if (norm <= 0.0) {
      // No neighbors in kernel support — fall back to nearest particle.
      int best = 0; double best_d2 = 1e30;
      for (size_t k = 0; k < p.size(); ++k) {
        double dx = p[k].r(0) - xj, dy = p[k].r(1) - yj;
        double dz = (D == 3) ? (p[k].r(2) - etj) : 0.0;
        double d2 = dx*dx + dy*dy + dz*dz;
        if (d2 < best_d2) { best_d2 = d2; best = (int)k; }
      }
      return std::make_pair(p[best].p() / p[best].T(),
                            p[best].T() * constants::hbarc_MeVfm);
    }
    return std::make_pair(rho_sum / norm, T_sum / norm);
  };

  const std::string& mode = settingsPtr->jets_input_mode;
  if (mode == "file") {
    std::ifstream fin(settingsPtr->jets_input_file);
    if (!fin.is_open()) {
      std::cerr << "BBMG: cannot open jets_input_file '"
                << settingsPtr->jets_input_file << "'\n";
      std::abort();
    }
    std::string line; int pid = 0;
    while (std::getline(fin, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream iss(line);
      double xj, yj, etj, pTj, phij, rapidity_;
      if (!(iss >> xj >> yj >> etj >> pTj >> phij >> rapidity_)) continue;
      auto [rho0, T_MeV] = sph_average(xj, yj, etj);
      push_pair(xj, yj, etj, phij, rapidity_, rho0, T_MeV, pid++);
    }
  }
  else if (mode == "hardcoded") {
    auto [rho0, T_MeV] = sph_average(settingsPtr->jets_x0,
                                     settingsPtr->jets_y0,
                                     settingsPtr->jets_eta0);
    push_pair(settingsPtr->jets_x0, settingsPtr->jets_y0, settingsPtr->jets_eta0,
              settingsPtr->jets_phi, settingsPtr->jets_rapidity,
              rho0, T_MeV, 0);
  }
  else if (mode == "sampled") {
    // Sample njet jets: position from a hot SPH particle, phi uniform in [0, 2π),
    // rapidity uniform in [-rapidity_max, +rapidity_max] (or fixed if sample_rapidity=false).
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> uphi(0.0, 2.0 * PI);
    std::uniform_real_distribution<double> uy(-settingsPtr->jets_rapidity_max,
                                               settingsPtr->jets_rapidity_max);
    std::uniform_int_distribution<>        uisph(0, (int)p.size() - 1);
    for (int i = 0; i < settingsPtr->jets_njet; ++i) {
      int rsp;
      do { rsp = uisph(gen); }
      while (p[rsp].T() * constants::hbarc_MeVfm < Freezeout_Temp);
      double y_s = settingsPtr->jets_sample_rapidity ? uy(gen) : settingsPtr->jets_rapidity;
      push_pair(p[rsp].r(0), p[rsp].r(1),
                (D == 3) ? p[rsp].r(2) : 0.0,
                uphi(gen), y_s,
                p[rsp].p() / p[rsp].T(),
                p[rsp].T() * constants::hbarc_MeVfm, i);
    }
  }
  else {
    std::cerr << "BBMG: unknown jets_input_mode '" << mode
              << "' (expected 'file', 'hardcoded', or 'sampled')\n";
    std::abort();
  }

  num_jets = (int)jetInfo_host.size();
  jetInfo = Cabana::AoSoA<jets, DeviceType, VECTOR_LENGTH>("jetInfo", num_jets);
  jetFreezeOut = jetInfo_host;
  copy_host_to_device_BBMG();
}
//////////////////////////////////////////////////////


//ccake::BBMG::flow(field &f) { return f.gam*(1-f.vmag*cos(f.phi-f.vang)); }

/*double BBMG::efluc()
{
  int random_variable = std::rand()/RAND_MAX;
  double zeta         = random_variable*(q+2.);
  return (1.+q) / pow(q+2, 1+q) * pow(q+2.-zeta, q);
}*/
template <unsigned int D>
void ccake::BBMG<D>::copy_host_to_device_BBMG(){

  using SerialHost = Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>;
  const int N = (int)jetInfo_host.size();
  Cabana::AoSoA<ccake::jets, SerialHost, 8> jets_h("jets_h", N);
  jets_VIEW(host_, jets_h);

  for (int ijet=0; ijet < N; ++ijet){
    for (int i=0; i<(int)D; ++i)
        host_position(ijet, i) = jetInfo_host[ijet].r[i];
    host_rho0(ijet) = jetInfo_host[ijet].rho0;
    host_x(ijet) = jetInfo_host[ijet].x;
    host_y(ijet) = jetInfo_host[ijet].y;
    host_eta(ijet) = jetInfo_host[ijet].eta;
    host_rho(ijet) = jetInfo_host[ijet].rho;
    host_T(ijet) = jetInfo_host[ijet].T;
    host_T0(ijet) = jetInfo_host[ijet].T0;
    host_phi(ijet) = jetInfo_host[ijet].phi;
    host_rapidity(ijet) = jetInfo_host[ijet].rapidity;
    host_energy_lost(ijet) = jetInfo_host[ijet].energy_lost;
    host_e_loss(ijet) = jetInfo_host[ijet].e_loss;
    host_gam(ijet) = jetInfo_host[ijet].gam;
    host_vmag(ijet) = jetInfo_host[ijet].vmag;
    host_vang(ijet) = jetInfo_host[ijet].vang;
    host_PID(ijet) = jetInfo_host[ijet].PID;
    host_Frozen(ijet) = jetInfo_host[ijet].Frozen;
  }

  Cabana::deep_copy( jetInfo, jets_h );
  Kokkos::fence();

}

template <unsigned int D>
void BBMG<D>::copy_device_to_host_BBMG()
{


  #ifdef DEBUG
  formatted_output::detail("Copying data from device to host");
  #endif
  //Auxiliary AoSoA for copying particles from/to host
  auto jets_h = Cabana::create_mirror_view_and_copy(Kokkos::HostSpace(), jetInfo);
  Kokkos::fence();
  jets_VIEW(host_, jets_h);
  const int N = (int)jetInfo_host.size();
  for (int ijet=0; ijet < N; ++ijet){
    for (int i=0; i<(int)D; ++i) jetInfo_host[ijet].r[i] = host_position(ijet, i);
    jetInfo_host[ijet].x = host_x(ijet);
    jetInfo_host[ijet].y = host_y(ijet);
    jetInfo_host[ijet].eta = host_eta(ijet);
    jetInfo_host[ijet].rho0 = host_rho0(ijet);
    jetInfo_host[ijet].rho  = host_rho(ijet);
    jetInfo_host[ijet].T = host_T(ijet);
    jetInfo_host[ijet].phi  = host_phi(ijet);
    jetInfo_host[ijet].rapidity = host_rapidity(ijet);
    jetInfo_host[ijet].energy_lost = host_energy_lost(ijet);
    jetInfo_host[ijet].e_loss = host_e_loss(ijet);
    jetInfo_host[ijet].gam  = host_gam(ijet);
    jetInfo_host[ijet].vmag = host_vmag(ijet);
    jetInfo_host[ijet].vang  = host_vang(ijet);
    jetInfo_host[ijet].PID = host_PID(ijet);
  }
}

template <unsigned int D>
void BBMG<D>::propagate()
{
  double tau  = systemPtr->t;

  jets_VIEW(jet_, jetInfo); 


  //#pragma omp parallel for schedule(auto)
  for (int i=0; i < jetInfo.size(); i++)  ///\@todo: Loop over jet_array index
  {
    jet_e_loss(i) = 0.0; // reset; set below only for active jets
    if (jet_T(i) > 150)
    {
    // Propagation direction for a massless jet with rapidity y in Milne.
    // With Δy ≡ y − η_s:  dx/dτ = cos(φ)/cosh(Δy), dy/dτ = sin(φ)/cosh(Δy),
    //                     dη_s/dτ = tanh(Δy)/τ.
    // In 2D, Δy = 0 → cosh=1, tanh=0, recovering pure transverse motion.
    const double dy   = (D == 3) ? jet_rapidity(i) - jet_position(i, 2) : 0.0;
    const double sech = 1.0 / cosh(dy);
    const double dir[3] = { cos(jet_phi(i)) * sech,
                            sin(jet_phi(i)) * sech,
                            (D == 3) ? tanh(dy) / tau : 0.0 };
    for (int idir = 0; idir < (int)D; ++idir)
      jet_position(i, idir) += vjet * settingsPtr->dt * dir[idir];

    Kokkos::View<double*> norm("norm",1);
    Kokkos::deep_copy(norm,0.0);
    double den = 0, den2 = 0;
    jet_T(i) = 0;
    jet_rho(i) = 0;
    for (int idir = 0; idir < (int)D; ++idir)
      jet_velocity(i, idir) = 0;

    CREATE_VIEW(device_, systemPtr->cabana_particles);
    auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());
    auto jet_interpolation = KOKKOS_LAMBDA(const int is, const int ia){
    // Always-factorized for D=3; falls back to SPHkernel<D> for D=1,2.
    double r_sph[D], r_jet[D];
    for (int idir = 0; idir < (int)D; ++idir) {
      r_sph[idir] = device_position.access(is, ia, idir);
      r_jet[idir] = jet_position(i, idir);
    }
    double kern = SPHkernel<D>::kernel(r_sph, r_jet, settingsPtr->hT, settingsPtr->hT);
    if (kern > 0.0)
    {
      Kokkos::atomic_add( &norm(0), kern);
      Kokkos::atomic_add( &jet_T(i), device_thermo.access(is, ia, thermo_info::T) * constants::hbarc_MeVfm * kern);
      Kokkos::atomic_add( &jet_rho(i), (device_thermo.access(is, ia, thermo_info::p)/device_thermo.access(is, ia, thermo_info::T)) * kern);

      for (int idir = 0; idir < (int)D; ++idir)
        Kokkos::atomic_add(&jet_velocity(i, idir), device_hydro_vector.access(is, ia, ccake::hydro_info::v, idir) * kern);

      }
   };
    ///////////////////////////////////////////////////////////////////
    Cabana::simd_parallel_for(simd_policy, jet_interpolation, "jet_interpolation");
    jet_T(i) /= norm(0);
    jet_rho(i) /= norm(0);
    for (int idir = 0; idir < (int)D; ++idir)
      jet_velocity(i, idir) /= norm(0);

    // |v_fluid|² with Milne metric (η component carries τ² in 3D hyperbolic).
    // make_covariant handles the D=2 vs D=3 sign/scaling internally.
    double t2 = 1.0;
    if constexpr (D == 3) {
      if (settingsPtr->coordinate_system == "hyperbolic") t2 = tau*tau;
    }
    milne::Vector<double, D> v_fluid;
    for (int idir = 0; idir < (int)D; ++idir) v_fluid(idir) = jet_velocity(i, idir);
    milne::Vector<double, D> v_fluid_cov = v_fluid;
    v_fluid_cov.make_covariant(t2);
    double vmag_sq = -milne::contract(v_fluid, v_fluid_cov);
    jet_vmag(i)    = sqrt(vmag_sq);
    if (jet_vmag(i) > 1) {
      cout << "The magnitude of the velocity returned greater than c, something is wrong. ABORTING" << endl;
      abort();
    }
    jet_vang(i) = atan2(jet_velocity(i, 1), jet_velocity(i, 0));
    jet_gam(i)  = 1.0 / sqrt(1.0 - vmag_sq);

    // Flow factor: γ_fluid · (1 − v⃗_jet · v⃗_fluid). η component carries τ².
    double v_dot = jet_velocity(i, 0) * dir[0]
                 + jet_velocity(i, 1) * dir[1];
    if constexpr (D == 3)
      v_dot += t2 * jet_velocity(i, 2) * dir[2];
    jet_flow(i) = jet_gam(i) * (1.0 - v_dot);
    double kappa = get_kappa(jet_T(i) / 1000); //The /1000 here is to move temps from MeV to GeV to follow Barbara's plot, same as above

    // Per-step energy loss ΔE = κ·τ^z·ρ^c·dt·flow; accumulate in the line
    // integral and store as e_loss for the source-term back-reaction.
    double dL =  kappa * exp(z*log(tau)) * exp(c*log(jet_rho(i))) * settingsPtr->dt * jet_flow(i);
    Kokkos::atomic_add( &jet_energy_lost(i), dL );
    jet_e_loss(i) = dL;

    }
  }
}

 //ccake namespace
#endif


