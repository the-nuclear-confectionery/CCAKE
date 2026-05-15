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
  using jets = Cabana::MemberTypes<double[2], // position
                                  double[2], //velocity
                                  double, // initial x
                                  double, // initial y
                                  double, // rho0
                                  double,    // rho
                                  double,    // Temp
                                  double,    // Initial Temp
                                  double,    // phi
                                  double,    // line integral
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
      rho0,
      rho,
      T,
      T0,
      phi,
      line_int,
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
  auto CONCAT(prefix, rho0) = Cabana::slice<jets_enum::rho0>(jets_aosoa); \
  auto CONCAT(prefix, rho) = Cabana::slice<jets_enum::rho>(jets_aosoa); \
  auto CONCAT(prefix, T) = Cabana::slice<jets_enum::T>(jets_aosoa); \
  auto CONCAT(prefix, T0) = Cabana::slice<jets_enum::T0>(jets_aosoa); \
  auto CONCAT(prefix, phi) = Cabana::slice<jets_enum::phi>(jets_aosoa); \
  auto CONCAT(prefix, line_int) = Cabana::slice<jets_enum::line_int>(jets_aosoa); \
  auto CONCAT(prefix, gam) = Cabana::slice<jets_enum::gam>(jets_aosoa); \
  auto CONCAT(prefix, vmag) = Cabana::slice<jets_enum::vmag>(jets_aosoa); \
  auto CONCAT(prefix, vang) = Cabana::slice<jets_enum::vang>(jets_aosoa); \
  auto CONCAT(prefix, flow) = Cabana::slice<jets_enum::flow>(jets_aosoa); \
  auto CONCAT(prefix, PID) = Cabana::slice<jets_enum::PID>(jets_aosoa); \
  auto CONCAT(prefix, Frozen) = Cabana::slice<jets_enum::Frozen>(jets_aosoa);


template <unsigned int D>
// I HAVE WRITTEN A 2 IN ANY PLACE I WOULD DEFINE D DIMENSIONS
// FIX IT EVENTUALLY WHEN BBMG DEPENDS ON D
// THERE IS ONE IN OUTPUT.H, FIX IT
// THERE IS ONE IN OUTPUT.CPP, FIX IT

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

    int phi_bins = systemPtr->jets_phi_bins;
    const static int phimax = 7;
    double phi[phimax];
    

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
    void initial_one_jet();

    Cabana::AoSoA<jets, DeviceType,VECTOR_LENGTH> jetInfo;

    struct field  // All fetched info for jet initialization and propagation
    {
        int sph, on;
        double rho, rho0, T, v[2];
        double r[2], phi, line_int, x, y;
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

template <unsigned int D>
double BBMG<D>::get_kappa(double T) //This is taken from one plot we found of Barbara's, no reason to believe this should be set in stone. Will test multiple Kappa's
{ 
  return 2.5*exp(-9.22*T) + 0.04;
}

template <unsigned int D>
void BBMG<D>::initial()
{
      srand( time( NULL ) );
  Freezeout_Temp    = 150;// Decoherence temperature in MeV
  Cg                = 3; // Cassimir const gluons
  Cq                = 4./3; // Cassimir const quarks
  // All next quantities are part of the BBMG parameters
  q                 = systemPtr->jets_Fluctuations; // Fluctuation parameter, subject to change
  z                 = systemPtr->jets_Length_scaling; // path length dependence
  a                 = systemPtr->jets_Energy_scaling; //Initial jet energy dependence
  c                 = (2+z-a)/3; //medium temperature dependence
  num_jets          = 200000; // Number of jets per event for oversampling
  //phimax            = 14;
  
  //===============================================
  vjet  = 1;

  gridx = settingsPtr->stepx;
  gridy = settingsPtr->stepy;

  for (int i = 0; i < phi_bins; i++)
  {
    phi[i] = i*PI/phi_bins;
  }
  //Setting final energy as a start point for the integration; This is starting in GeV
  Pfg = 10;
  Pfq = 10;

  Pfg /= constants::hbarc_GeVfm; //Converting into femtometers
  Pfq /= constants::hbarc_GeVfm;

////////////////////////////////////////////////////////////////////////////////////////////

    //vector<Particle>  p_bbmg;

    jetInfo = Cabana::AoSoA<jets, DeviceType,VECTOR_LENGTH>("jetInfo", 2*num_jets);
    auto & p = systemPtr->particles;

    int back_to_back = 2;
    std::random_device rd; // For true randomness
    std::mt19937 gen(rd()); // Mersenne Twister generator
    std::uniform_int_distribution<> dis(0, p.size() - 1); // Distribute over the valid range
    int random_sph_particle;
    for (int i = 0; i < num_jets; ++i)
    {
        // While loop to check temperature of chosen sph particle with if statement
        
        bool above_freeze = true;
        
        while (above_freeze)
          {
            random_sph_particle = dis(gen);
            double Temp = p[random_sph_particle].T();
            if (Temp * constants::hbarc_MeVfm >= Freezeout_Temp)
              above_freeze = false;
            else
              above_freeze = true;
          }
        
        
        
        //Density from pressure over temperature
        double rsub = p[random_sph_particle].p() / p[random_sph_particle].T();
        field sph_particle; //field of all sph particles where we take necessary line integral info
        sph_particle.r[0] = p[random_sph_particle].r(0);
        sph_particle.r[1] = p[random_sph_particle].r(1);
        sph_particle.x    = p[random_sph_particle].r(0);
        sph_particle.y    = p[random_sph_particle].r(1);
        sph_particle.rho0 = rsub; //Density left in terms of femtometers
        sph_particle.T = p[random_sph_particle].T() * constants::hbarc_MeVfm;
        //sph_particle.T0 = p[random_sph_particle].T() * constants::hbarc_MeVfm;
        //The above line is to create histograms to compare with initial temperature distribution
        
        double kappa = get_kappa(sph_particle.T / 1000);

        sph_particle.line_int = 0.5 * kappa * exp(z * log(settingsPtr->t0)) * exp(c * log(sph_particle.rho0)) * settingsPtr->dt; // only if initial flow=0
        int phidist = rand() % phi_bins;
        for (int j = 0; j < back_to_back; j++)
        {
            sph_particle.phi = phi[phidist] + j * PI;
            sph_particle.PID = phidist + phi_bins * j;
            sph_particle.Frozen = 1;
            jetInfo_host.push_back(sph_particle);
        }
    }
    jetFreezeOut = jetInfo_host;


    copy_host_to_device_BBMG();
}
//////////////////////////////////////////////////////

template <unsigned int D>
void BBMG<D>::initial_one_jet()
{
      srand( time( NULL ) );
  Freezeout_Temp    = 150;// Decoherence temperature in MeV
  Cg                = 3; // Cassimir const gluons
  Cq                = 4./3; // Cassimir const quarks
  // All next quantities are part of the BBMG parameters
  q                 = 0; // Fluctuation parameter, subject to change
  z                 = 1; // path length dependence
  a                 = 0; //Initial jet energy dependence
  c                 = (2+z-a)/3; //medium temperature dependence
  num_jets          = 1;
  //phimax            = 14;
  
  //===============================================
  vjet  = 1;

  gridx = settingsPtr->stepx;
  gridy = settingsPtr->stepy;

  for (int i = 0; i < phi_bins; i++)
  {
    phi[i] = i*PI/phi_bins;
  }
  //Setting final energy as a start point for the integration; This is starting in GeV
  Pfg = 10;
  Pfq = 10;

  Pfg /= constants::hbarc_GeVfm; //Converting into femtometers
  Pfq /= constants::hbarc_GeVfm;

////////////////////////////////////////////////////////////////////////////////////////////

  jetInfo = Cabana::AoSoA<jets, DeviceType,VECTOR_LENGTH>("jetInfo", 2*num_jets);

    
  for ( int i = 0; i < systemPtr->particles.size(); ++i )
  {
    auto & p = systemPtr->particles[i];
    if (p.r(0)==1.71 && p.r(1)==0.51 && p.T()*constants::hbarc_MeVfm>=150)
      { 
        int index = i;
        cout << "sph particle index is " << index << endl;
        cout << "initial sph temperature is " << p.T()*constants::hbarc_MeVfm << endl;
        int back_to_back = 2;
        double rsub = p.p() / p.T();
        field sph_particle; //field of all sph particles where we take necessary line integral info
        sph_particle.r[0] = p.r(0);
        sph_particle.x = p.r(0);
        //cout << "Positions (x) of each particle in the grid is " << p.r(0) << "\n";
        sph_particle.r[1] = p.r(1);
        sph_particle.y = p.r(1);
        sph_particle.rho0 = rsub; //Density left in terms of femtometers
        //sph_particle.sph = i;
        sph_particle.T = p.T() * constants::hbarc_MeVfm;
        sph_particle.T0 = p.T() * constants::hbarc_MeVfm;
        //The above line is to create histograms to compare with initial temperature distribution

        double kappa = get_kappa(sph_particle.T / 1000);

        sph_particle.line_int = 0.5 * kappa * exp(z * log(settingsPtr->t0)) * exp(c * log(sph_particle.rho0)) * settingsPtr->dt; // only if initial flow=0
        for (int j = 0; j < back_to_back; j++)
        {
            sph_particle.phi = PI/7 + j * PI;
            sph_particle.PID = 1 + phi_bins * j;
            jetInfo_host.push_back(sph_particle);
        }
      }

    }
  copy_host_to_device_BBMG();
    //copy_device_to_host_BBMG();
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
  //Auxiliary AoSoA for copying particles from/to host
  Cabana::AoSoA<ccake::jets, SerialHost, 8> jets_h("jets_h", 2*num_jets);
  
  jets_VIEW(host_, jets_h);


  //Fill host arrays
  for (int ijet=0; ijet < 2*num_jets; ++ijet){
    //Copy hydro space matrix
    for (int i=0; i<2; ++i){
        host_position(ijet, i) = jetInfo_host[ijet].r[i];
          }
    host_rho0(ijet) = jetInfo_host[ijet].rho0;
    host_x(ijet) = jetInfo_host[ijet].x;
    host_y(ijet) = jetInfo_host[ijet].y;
    host_rho(ijet) = jetInfo_host[ijet].rho;
    host_T(ijet) = jetInfo_host[ijet].T;
    host_T0(ijet) = jetInfo_host[ijet].T0;
    host_phi(ijet) = jetInfo_host[ijet].phi;
    host_line_int(ijet) = jetInfo_host[ijet].line_int;
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
  for (int ijet=0; ijet < 2*num_jets; ++ijet){

    for (int i=0; i<2; ++i)
    {
      jetInfo_host[ijet].r[i] = host_position(ijet, i);
    }
    jetInfo_host[ijet].x = host_x(ijet);
    jetInfo_host[ijet].y = host_y(ijet);
    jetInfo_host[ijet].rho0 = host_rho0(ijet);
    jetInfo_host[ijet].rho  = host_rho(ijet);
    jetInfo_host[ijet].T = host_T(ijet);
    jetInfo_host[ijet].phi  = host_phi(ijet);
    jetInfo_host[ijet].line_int = host_line_int(ijet);
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
    if (jet_T(i) > 150)
    {
    // propagate x,y position of jet on top of sph particles
    jet_position(i,0) += vjet * settingsPtr->dt * cos(jet_phi(i)); //Flow is not here to alter the direction of the jet
    jet_position(i,1) += vjet * settingsPtr->dt * sin(jet_phi(i));

    //double flow;
    Kokkos::View<double*> norm("norm",1);
    Kokkos::deep_copy(norm,0.0);
    double den = 0, den2 = 0;
    jet_T(i) = 0;
    jet_rho(i) = 0;
    jet_velocity(i, 0) = 0;
    jet_velocity(i, 1) = 0;

    CREATE_VIEW(device_, systemPtr->cabana_particles);
    auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());
    auto jet_interpolation = KOKKOS_LAMBDA(const int is, const int ia){
    //double r[D];
    double dx    = device_position.access(is, ia, 0) - jet_position(i, 0);
    double dy    = device_position.access(is, ia, 1) - jet_position(i, 1);

    double rdiff = sqrt(dx*dx+dy*dy);
    if (rdiff<2*settingsPtr->hT)
    {
      double kern = SPHkernel<D>::kernel(rdiff,settingsPtr->hT);
      Kokkos::atomic_add( &norm(0), kern);
      Kokkos::atomic_add( &jet_T(i), device_thermo.access(is, ia, thermo_info::T) * constants::hbarc_MeVfm * kern);
      Kokkos::atomic_add( &jet_rho(i), (device_thermo.access(is, ia, thermo_info::p)/device_thermo.access(is, ia, thermo_info::T)) * kern);
    
      for (int idir = 0; idir < 2; ++idir)
        Kokkos::atomic_add(&jet_velocity(i, idir), device_hydro_vector.access(is, ia, ccake::hydro_info::v, idir) * kern);
  
      }
   };
    ///////////////////////////////////////////////////////////////////
    Cabana::simd_parallel_for(simd_policy, jet_interpolation, "jet_interpolation"); 
      // Put simd_parallel above this line
    jet_T(i) /= norm(0);
    jet_rho(i) /= norm(0);
    for (int idir = 0; idir < 2; ++idir)
      jet_velocity(i, idir) /= norm(0);

    jet_vmag(i) = sqrt(jet_velocity(i, 0) * jet_velocity(i, 0) + jet_velocity(i, 1) * jet_velocity(i, 1));
    if (jet_vmag(i)>1)
      {
        cout << "The magnitude of the velocity returned greater than c, something is wrong. ABORTING" << endl;
        abort();
      }
      jet_vang(i) = atan2( jet_velocity(i, 1), jet_velocity(i, 0));
      jet_gam(i) = 1.0 / sqrt( 1.0 - jet_vmag(i) * jet_vmag(i));


    jet_flow(i) = jet_gam(i)*(1-jet_vmag(i)*cos(jet_phi(i)-jet_vang(i)));
    double kappa = get_kappa(jet_T(i) / 1000); //The /1000 here is to move temps from MeV to GeV to follow Barbara's plot, same as above

    Kokkos::atomic_add( &jet_line_int(i), kappa * exp(z*log(tau)) * exp(c*log(jet_rho(i))) * settingsPtr->dt * jet_flow(i));
      
    }
  }
}

 //ccake namespace
#endif


