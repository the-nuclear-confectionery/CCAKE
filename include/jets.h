#ifndef BBMG_H
#define BBMG_H

//#include <cerrno>
//#include <cstdio>
//#include <cstdlib>
//#include <cstring>
//#include <ctime>
//#include <fstream>
//#include <iostream>
//#include <sstream>
//#include <string>
//#include <vector>
//#include <cmath>

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
//#include "table_delaunay.h"
#include "thermodynamic_info.h"
#include "milne.hpp"

#include <Cabana_Core.hpp>



//using namespace std;

namespace ccake
{
  using jets = Cabana::MemberTypes<double[2], // position
                                  double[2], //velocity
                                  double, // rho0
                                  double,    // rho
                                  double,    // Temp
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
      rho0,
      rho,
      T,
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
  //using jet_freeze = Cabana::AoSoA<jets, DeviceType,VECTOR_LENGTH>;
  
  #define jets_VIEW(prefix, jets_aosoa) \
  auto CONCAT(prefix, position) = Cabana::slice<jets_enum::position>(jets_aosoa); \
  auto CONCAT(prefix, velocity) = Cabana::slice<jets_enum::velocity>(jets_aosoa); \
  auto CONCAT(prefix, rho0) = Cabana::slice<jets_enum::rho0>(jets_aosoa); \
  auto CONCAT(prefix, rho) = Cabana::slice<jets_enum::rho>(jets_aosoa); \
  auto CONCAT(prefix, T) = Cabana::slice<jets_enum::T>(jets_aosoa); \
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

//friend class ccake::Output;

private:
    std::shared_ptr<Settings> settingsPtr;
    std::shared_ptr<SystemState<D>> systemPtr;

    int z, a, c; ///\@todo: Create an input for these parameters
    int num_jets;
    double Freezeout_Temp;
    double area;
    double vjet; // Taken to be c for jets
    double Cg, Cq, q;
    double gridx, gridy;
    
    double Pfg, Pfq; ///\@todo: Move or eliminate based on need for energy loss equation with energy dependence


    const static int phimax = 7;
    double phi[phimax];
    //double Rjetq[phimax], Rjetg[phimax];
    //vector<double> rr;
    
    //vector<double> RAAq;
    //vector<double> RAAg; 

    //jet_array jetInfo;
    //jet_array jetFreezeOut;
    

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


    // Testing to see if this being public allows deep copy to work and allows copy from device to host
    Cabana::AoSoA<jets, DeviceType,VECTOR_LENGTH> jetInfo;

    struct field  // All fetched info for jet initialization and propagation
    {
        int sph, on;
        double rho, rho0, T, v[2];
        double r[2], phi, line_int;
        int PID, Frozen;
        double gam, vmag, vang, flow;
        //double T0;
    };

    vector<field> jetInfo_host; ///\@todo: Change these vectors to the AoSoA's
    vector<field> jetFreezeOut;

    //jet_array jets;
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
  q                 = 0; // Fluctuation parameter, subject to change
  z                 = 1; // path length dependence
  a                 = 0; //Initial jet energy dependence
  c                 = (2+z-a)/3; //medium temperature dependence
  num_jets          = 5000; // Number of jets per event for oversampling
  //phimax            = 14;
  
  //===============================================
  vjet  = 1;
  // area is taken from parameter h read in from input parameters file
  //area  = PI*pow(2.*systemPtr->h,2);

  gridx = settingsPtr->stepx;
  gridy = settingsPtr->stepy;
  //cout << "Gridx and gridy are " << gridx << "," << gridy << endl << endl;

  for (int i = 0; i < phimax; i++)
  {
    phi[i] = i*PI/7;
  }
  //Setting final energy as a start point for the integration; This is starting in GeV
  Pfg = 10;
  Pfq = 10;

  Pfg /= constants::hbarc_GeVfm; //Converting into femtometers
  Pfq /= constants::hbarc_GeVfm;

////////////////////////////////////////////////////////////////////////////////////////////

    //vector<Particle>  p_bbmg;

    jetInfo = Cabana::AoSoA<jets, DeviceType,VECTOR_LENGTH>("jetInfo", 2*num_jets);
    //jetFreezeOut = jet_freeze("jetFreezeOut", 2*num_jets);
    //jetFreezeOut = jet_array("jetFreezeOut", 2*num_jets);
    auto & p = systemPtr->particles;
    
    /*for ( int i = 0; i < systemPtr->particles.size(); ++i )
    {
      auto & p = systemPtr->particles[i];
      if (p.T() * constants::hbarc_MeVfm >= Freezeout_Temp)
      {
        p_bbmg.push_back(p);
      }
    }*/

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
        sph_particle.rho0 = rsub; //Density left in terms of femtometers
        sph_particle.T = p[random_sph_particle].T() * constants::hbarc_MeVfm;
        //sph_particle.T0 = p[random_sph_particle].T() * constants::hbarc_MeVfm;
        //The above line is to create histograms to compare with initial temperature distribution
        
        double kappa = get_kappa(sph_particle.T / 1000);

        sph_particle.line_int = 0.5 * kappa * exp(z * log(settingsPtr->t0)) * exp(c * log(sph_particle.rho0)) * settingsPtr->dt; // only if initial flow=0
        int phidist = rand() % phimax;
        for (int j = 0; j < back_to_back; j++)
        {
            sph_particle.phi = phi[phidist] + j * PI;
            sph_particle.PID = phidist + phimax * j;
            sph_particle.Frozen = 1;
            jetInfo_host.push_back(sph_particle);
        }
    }
    jetFreezeOut = jetInfo_host;


    copy_host_to_device_BBMG();
    //copy_device_to_host_BBMG();
}


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
    host_rho(ijet) = jetInfo_host[ijet].rho;
    host_T(ijet) = jetInfo_host[ijet].T;
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
  cout << "jets_VIEW instantiation..." << endl;
  jets_VIEW(host_, jets_h);
  for (int ijet=0; ijet < 2*num_jets; ++ijet){
    //int id = host_id(ijet);
    for (int i=0; i<2; ++i)
    {
      jetFreezeOut[ijet].r[i] = host_position(ijet, i);
    }
    jetFreezeOut[ijet].rho0 = host_rho0(ijet);
    jetFreezeOut[ijet].rho  = host_rho(ijet);
    jetFreezeOut[ijet].T = host_T(ijet);
    jetFreezeOut[ijet].phi  = host_phi(ijet);
    jetFreezeOut[ijet].line_int = host_line_int(ijet);
    jetFreezeOut[ijet].gam  = host_gam(ijet);
    jetFreezeOut[ijet].vmag = host_vmag(ijet);
    jetFreezeOut[ijet].vang  = host_vang(ijet);
    jetFreezeOut[ijet].PID = host_PID(ijet);
}
}  

template <unsigned int D>
void BBMG<D>::propagate()
{
  double tau  = systemPtr->t;// + settingsPtr->t0;

  jets_VIEW(jet_, jetInfo);
  //jets_VIEW(jetFreeze_, jetFreeze);     


  //#pragma omp parallel for schedule(auto)
  for (int i=0; i < jetInfo.size(); i++)  ///\@todo: Loop over jet_array index 
  {
    if (jet_T(i) > 150)
    {
    // propagate x,y position of jet on top of sph particles
    jet_position(i,0) += vjet * settingsPtr->dt * cos(jet_phi(i)); //Flow is not here to alter the direction of the jet
    jet_position(i,1) += vjet * settingsPtr->dt * sin(jet_phi(i));
    //jetPropagation.r[0] += vjet * settingsPtr->dt * cos(jetPropagation.phi); //Flow is not here to alter the direction of the jet
    //jetPropagation.r[1] += vjet * settingsPtr->dt * sin(jetPropagation.phi);

    //double flow;
    Kokkos::View<double*> norm("norm",1);
    Kokkos::deep_copy(norm,0.0);
    //double norm = 0;
    double den = 0, den2 = 0;
    jet_T(i) = 0;
    jet_rho(i) = 0;
    jet_velocity(i, 0) = 0;
    jet_velocity(i, 1) = 0;
    //f.T    = 0;
    //f.rho  = 0;
    //f.v[0] = 0;
    //f.v[1] = 0;

    CREATE_VIEW(device_, systemPtr->cabana_particles);
    auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());
    auto jet_interpolation = KOKKOS_LAMBDA(const int is, const int ia){
    //double r[D];
    double dx    = device_position.access(is, ia, 0) - jet_position(i, 0);
    double dy    = device_position.access(is, ia, 1) - jet_position(i, 1);

    double rdiff = sqrt(dx*dx+dy*dy);
    if (rdiff<2*settingsPtr->hT)
    {
      //den++;
      //den2     += p.norm_spec.s;
      double kern = SPHkernel<D>::kernel(rdiff,settingsPtr->hT);
      Kokkos::atomic_add( &norm(0), kern);
      //norm         += kern;
      Kokkos::atomic_add( &jet_T(i), device_thermo.access(is, ia, thermo_info::T) * constants::hbarc_MeVfm * kern);
      //jet_T(i) += device_thermo.access(is, ia, thermo_info::T) * ccake::constants::hbarc_MeVfm * kern;
      //f.T      += p.T()*ccake::constants::hbarc_MeVfm*kern;
      //cout << "Interpolated temp value is " << f.T << "\n";
      Kokkos::atomic_add( &jet_rho(i), (device_thermo.access(is, ia, thermo_info::p)/device_thermo.access(is, ia, thermo_info::T)) * kern);
      
      //jet_rho(i) += (device_thermo.access(is, ia, ccake::thermo_info::p)/device_thermo.access(is, ia, ccake::thermo_info::T)) * kern
      //f.rho    += (p.p()/p.T())*kern;
      for (int idir = 0; idir < 2; ++idir)
        Kokkos::atomic_add(&jet_velocity(i, idir), device_hydro_vector.access(is, ia, ccake::hydro_info::v, idir) * kern);
        //if (device_hydro_vector.access(is, ia, ccake::hydro_info::v, idir)* kern)
        //jet_velocity(i, idir) += device_hydro_vector.access(is, ia, ccake::hydro_info::v, idir);
      
      //f.v[0]   += p.hydro.v(0)*kern;
      //f.v[1]   += p.hydro.v(1)*kern;
           }
          };
    ///////////////////////////////////////////////////////////////////
    Cabana::simd_parallel_for(simd_policy, jet_interpolation, "jet_interpolation"); 
      // Put simd_parallel above this line
    jet_T(i) /= norm(0);
    // f.T       /= norm;
    jet_rho(i) /= norm(0);
    // f.rho     /= norm;
    for (int idir = 0; idir < 2; ++idir)
      jet_velocity(i, idir) /= norm(0);
    // f.v[0]    /= norm;
    // f.v[1]    /= norm;
    jet_vmag(i) = sqrt(jet_velocity(i, 0) * jet_velocity(i, 0) + jet_velocity(i, 1) * jet_velocity(i, 1));
    //cout << "This is the magnitude of the velocity around each jet: " << jet_vmag(i) << endl;
    //cout << "These are the components in x: " << jet_velocity(i, 0) << endl;
    //cout << "These are the components in y: " << jet_velocity(i, 1) << endl;
    //f.vmag     = sqrt( f.v[0]*f.v[0] + f.v[1]*f.v[1] );
    if (jet_vmag(i)>1)
      {
        cout << "The magnitude of the velocity returned greater than c, something is wrong. ABORTING" << endl;
        abort();
      }
      jet_vang(i) = atan2( jet_velocity(i, 1), jet_velocity(i, 0));
      //f.vang     = atan2( f.v[1], f.v[0] );
      jet_gam(i) = 1.0 / sqrt( 1.0 - jet_vmag(i) * jet_vmag(i));
      //f.gam      = 1.0 / sqrt( 1.0 - f.vmag*f.vmag );
    }

    
    
    jet_flow(i) = jet_gam(i)*(1-jet_vmag(i)*cos(jet_phi(i)-jet_vang(i)));
    double kappa = get_kappa(jet_T(i) / 1000); //The /1000 here is to move temps from MeV to GeV to follow Barbara's plot, same as above

    Kokkos::atomic_add( &jet_line_int(i), kappa * exp(z*log(tau)) * exp(c*log(jet_rho(i))) * settingsPtr->dt * jet_flow(i));
    //jet_line_int(i) += kappa * exp(z*log(tau)) * exp(c*log(jet_rho(i))) * settingsPtr->dt * flow;
    //jetPropagation.line += kappa * exp(z*log(tau)) * exp(c*log(jetPropagation.rho)) * settingsPtr->dt * flow(jetPropagation);
    //countyes++;
    //cout << "Jet directions still going: " << jetPropagation.phi << endl;
    //inter( jetPropagation ); //interpolation of the area around the jet
  }
}
        /*auto condition = [this](auto& jetPropagation) { ///\@todo: Get rid of the removal routine as we have two arrays that contain both frozen and non-frozen jets
            return jetPropagation.T <= Freezeout_Temp;
        };

        auto new_end = std::remove_if(jetInfo.begin(), jetInfo.end(),
            [this, &condition](auto& jetPropagation) {
                if (condition(jetPropagation)) {
                    jetFreezeOut.push_back(jetPropagation);
                    return true; // Mark element for removal
                }
                return false; // Keep element
            });


        // Erase the removed elements from the source vector
        jetInfo.erase(new_end, jetInfo.end());*/




/*void BBMG::inter( field &f ) 
{
  double norm = 0;
  double den = 0, den2 = 0;
  f.T    = 0;
  f.rho  = 0;
  f.v[0] = 0;
  f.v[1] = 0;
  #pragma omp parallel for schedule(auto)
  for ( auto & p : systemPtr->particles )
  {
    double dx    = p.r(0)-f.r[0];
    double dy    = p.r(1)-f.r[1];

    double rdiff = sqrt(dx*dx+dy*dy)/systemPtr->h;
    

    if (rdiff<2)
    {
      den++;
      den2     += p.norm_spec.s;
      double kern = kernel::kernel(rdiff);
      norm         += kern;
      f.T      += p.T()*constants::hbarc_MeVfm*kern;
      //cout << "Interpolated temp value is " << f.T << "\n";
      f.rho    += (p.p()/p.T())*kern;
      f.v[0]   += p.hydro.v(0)*kern;
      f.v[1]   += p.hydro.v(1)*kern;
      
    }
  }
  f.T       /= norm;
  f.rho     /= norm;
  f.v[0]    /= norm;
  f.v[1]    /= norm;
  f.vmag     = sqrt( f.v[0]*f.v[0] + f.v[1]*f.v[1] );
  if (f.vmag>1)
  {
    cout << "The magnitude of the velocity returned greater than c, something is wrong. ABORTING" << endl;
    abort();
  }
  f.vang     = atan2( f.v[1], f.v[0] );
  f.gam      = 1.0 / sqrt( 1.0 - f.vmag*f.vmag );
}*/



 //ccake namespace
#endif

