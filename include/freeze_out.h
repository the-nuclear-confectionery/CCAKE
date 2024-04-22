#ifndef FREEZEOUT_H
#define FREEZEOUT_H

#include <algorithm>
#include <string>
#include <vector>

#include <Cabana_Core.hpp>

#include "eos.h"
#include "particle.h"
#include "system_state.h"
#include "constants.h"
#include "milne.hpp"


namespace ccake
{

  #ifndef CONCAT
    #define CONCAT_(a, b) a##b
    #define CONCAT(a, b) CONCAT_(a, b)
  #endif

  using FRZ = Cabana::MemberTypes<double[4][4], double[3], double[3], double[3], double[3], // shear, r, u, gradP, gradE
                                  double, double, double, // t, s, e
                                  double, double, double, // rhoB, rhoS, rhoQ
                                  double, double, double, double, // T, muB, muS, muQ
                                  double, double, double, // theta, bulk, sigma
                                  double, double, // shear33, inside
                                  double, double // wfz, cs2fz
                                  >; // 41*8 = 328 bytes per particle. 100K Particles = 32.8 MB
  using ResultsTypes = Cabana::MemberTypes<double[4][4], // shearsub
                                           double[3], double[3], double[3], //divT, divE, divP
                                           double[3], double[3], //rsub, uout
                                           double, double, double, // divTtemp, divPpress, divEener
                                           double, double, double, // gsub, bulksub, thetasub
                                           double, double, double, // swsub, shear33sub, tlist
                                           double, double, double, double, // sFO, Efluc, Tfluc, muBfluc
                                           double, double, double, // muSfluc, muQfluc, cs2fzfluc
                                           double, bool    // wfzfluc
                                           >; // 36*8 = 288 bytes per particle. 100K Particles = 28.8 MB
  namespace FRZ_enum{
    enum FRZ_members
    {
      shear, r, u, gradP, gradE,
      t, s, e,
      rhoB, rhoS, rhoQ,
      T, muB, muS, muQ,
      theta, bulk, sigma,
      shear33, inside,
      wfz, cs2fz
    };
  }
  namespace ResultsTypes_enum{
    enum ResultsTypes_members
    {
      shearsub,
      divT, divE, divP,
      rsub, uout,
      divTtemp, divPpress, divEener,
      gsub, bulksub, thetasub,
      swsub, shear33sub, tlist,
      sFO, Efluc, Tfluc, muBfluc,
      muSfluc, muQfluc, cs2fzfluc,
      wfzfluc, print
    };
  }

  using freeze_array = Cabana::AoSoA<FRZ,DeviceType,VECTOR_LENGTH>;
  using freeze_results = Cabana::AoSoA<ResultsTypes,DeviceType,VECTOR_LENGTH>;
  using mask_type = Cabana::AoSoA<Cabana::MemberTypes<bool>,DeviceType,VECTOR_LENGTH>;

  #define FRZ_VIEW(prefix,FRZ_aosoa) \
  auto CONCAT(prefix, shear) = Cabana::slice<FRZ_enum::shear>(FRZ_aosoa); \
  auto CONCAT(prefix, r) = Cabana::slice<FRZ_enum::r>(FRZ_aosoa); \
  auto CONCAT(prefix, u) = Cabana::slice<FRZ_enum::u>(FRZ_aosoa); \
  auto CONCAT(prefix, gradP) = Cabana::slice<FRZ_enum::gradP>(FRZ_aosoa); \
  auto CONCAT(prefix, gradE) = Cabana::slice<FRZ_enum::gradE>(FRZ_aosoa); \
  auto CONCAT(prefix, t) = Cabana::slice<FRZ_enum::t>(FRZ_aosoa); \
  auto CONCAT(prefix, s) = Cabana::slice<FRZ_enum::s>(FRZ_aosoa); \
  auto CONCAT(prefix, e) = Cabana::slice<FRZ_enum::e>(FRZ_aosoa); \
  auto CONCAT(prefix, rhoB) = Cabana::slice<FRZ_enum::rhoB>(FRZ_aosoa); \
  auto CONCAT(prefix, rhoS) = Cabana::slice<FRZ_enum::rhoS>(FRZ_aosoa); \
  auto CONCAT(prefix, rhoQ) = Cabana::slice<FRZ_enum::rhoQ>(FRZ_aosoa); \
  auto CONCAT(prefix, T) = Cabana::slice<FRZ_enum::T>(FRZ_aosoa); \
  auto CONCAT(prefix, muB) = Cabana::slice<FRZ_enum::muB>(FRZ_aosoa); \
  auto CONCAT(prefix, muS) = Cabana::slice<FRZ_enum::muS>(FRZ_aosoa); \
  auto CONCAT(prefix, muQ) = Cabana::slice<FRZ_enum::muQ>(FRZ_aosoa); \
  auto CONCAT(prefix, theta) = Cabana::slice<FRZ_enum::theta>(FRZ_aosoa); \
  auto CONCAT(prefix, bulk) = Cabana::slice<FRZ_enum::bulk>(FRZ_aosoa); \
  auto CONCAT(prefix, sigma) = Cabana::slice<FRZ_enum::sigma>(FRZ_aosoa); \
  auto CONCAT(prefix, shear33) = Cabana::slice<FRZ_enum::shear33>(FRZ_aosoa); \
  auto CONCAT(prefix, inside) = Cabana::slice<FRZ_enum::inside>(FRZ_aosoa); \
  auto CONCAT(prefix, wfz) = Cabana::slice<FRZ_enum::wfz>(FRZ_aosoa); \
  auto CONCAT(prefix, cs2fz) = Cabana::slice<FRZ_enum::cs2fz>(FRZ_aosoa);

  #define FRZ_RESULTS_VIEW(prefix,results_aosoa) \
  auto CONCAT(prefix, shearsub) = Cabana::slice<ResultsTypes_enum::shearsub>(results_aosoa); \
  auto CONCAT(prefix, divT) = Cabana::slice<ResultsTypes_enum::divT>(results_aosoa); \
  auto CONCAT(prefix, divE) = Cabana::slice<ResultsTypes_enum::divE>(results_aosoa); \
  auto CONCAT(prefix, divP) = Cabana::slice<ResultsTypes_enum::divP>(results_aosoa); \
  auto CONCAT(prefix, rsub) = Cabana::slice<ResultsTypes_enum::rsub>(results_aosoa); \
  auto CONCAT(prefix, uout) = Cabana::slice<ResultsTypes_enum::uout>(results_aosoa); \
  auto CONCAT(prefix, divTtemp) = Cabana::slice<ResultsTypes_enum::divTtemp>(results_aosoa); \
  auto CONCAT(prefix, divPpress) = Cabana::slice<ResultsTypes_enum::divPpress>(results_aosoa); \
  auto CONCAT(prefix, divEener) = Cabana::slice<ResultsTypes_enum::divEener>(results_aosoa); \
  auto CONCAT(prefix, gsub) = Cabana::slice<ResultsTypes_enum::gsub>(results_aosoa); \
  auto CONCAT(prefix, bulksub) = Cabana::slice<ResultsTypes_enum::bulksub>(results_aosoa); \
  auto CONCAT(prefix, thetasub) = Cabana::slice<ResultsTypes_enum::thetasub>(results_aosoa); \
  auto CONCAT(prefix, swsub) = Cabana::slice<ResultsTypes_enum::swsub>(results_aosoa); \
  auto CONCAT(prefix, shear33sub) = Cabana::slice<ResultsTypes_enum::shear33sub>(results_aosoa); \
  auto CONCAT(prefix, tlist) = Cabana::slice<ResultsTypes_enum::tlist>(results_aosoa); \
  auto CONCAT(prefix, sFO) = Cabana::slice<ResultsTypes_enum::sFO>(results_aosoa); \
  auto CONCAT(prefix, Efluc) = Cabana::slice<ResultsTypes_enum::Efluc>(results_aosoa); \
  auto CONCAT(prefix, Tfluc) = Cabana::slice<ResultsTypes_enum::Tfluc>(results_aosoa); \
  auto CONCAT(prefix, muBfluc) = Cabana::slice<ResultsTypes_enum::muBfluc>(results_aosoa); \
  auto CONCAT(prefix, muSfluc) = Cabana::slice<ResultsTypes_enum::muSfluc>(results_aosoa); \
  auto CONCAT(prefix, muQfluc) = Cabana::slice<ResultsTypes_enum::muQfluc>(results_aosoa); \
  auto CONCAT(prefix, cs2fzfluc) = Cabana::slice<ResultsTypes_enum::cs2fzfluc>(results_aosoa); \
  auto CONCAT(prefix, wfzfluc) = Cabana::slice<ResultsTypes_enum::wfzfluc>(results_aosoa); \
  auto CONCAT(prefix, print) = Cabana::slice<ResultsTypes_enum::print>(results_aosoa);


//TODO: This needs to be documented and split in header and
//source code

//==============================================================================
// A DICTIONARY FOR CONFUSINGLY NAMED QUANTITIES
//------------------------------------------------------------------------------
// tau    - current tau
// taup   - tau of previous timestep
// taupp  - tau of two timesteps ago
// Freeze - current freeze out status of particle (0 == freeze-out not begun,
//                                                 1 == freeze-out begun,
//                                                 3 == freeze-out basically done,
//                                                 4 == completely frozen out)
// frz1   - snapshot of particle state at taup
// frz2   - snapshot of particle state at taupp
// fback1 - snapshot of particle state at taup
// fback2 - snapshot of particle state at taupp
// fback3 - snapshot of particle state at tauppp(?)
// fback4 - snapshot of particle state at taupppp(?)
// cf     - basically the same as curfrz
// btrack - counts the number of nearest neighbors a particle has
//==============================================================================

template <unsigned int D>
class FreezeOut
{
  private:

    std::shared_ptr<Settings> settingsPtr    = nullptr;
    std::shared_ptr<SystemState<D>> systemPtr = nullptr;

    // freeze out structures
    freeze_array frz1;   ///< freeze out snapshots at last time step
    freeze_array frz2;   ///< freeze out snapshots at penultime time step
    freeze_array fback;  ///< Backup snapshots for cases of few neighbors
    freeze_array fback2; ///< Backup snapshots for cases of few neighbors
    freeze_array fback3; ///< Backup snapshots for cases of few neighbors
    freeze_array fback4; ///< Backup snapshots for cases of few neighbors
    mask_type mask;      ///< Mask for tagging particles which should not be copied from one snapshot to another

    //==========================================================================
    // freeze-out related quantities
    int cf            = 0;
    int frzc          = 0;
    double cs2        = 0;
    double tau        = 0; ///< Current time step
    double taup       = 0; ///< Previous time step
    double taupp      = 0; ///< Time step before previous
    double wfz        = 0; 


    vector<string> eosname;
    //==========================================================================

    /// @brief Initializes the freeze_out object with the given settings and system state.
    /// @details Access to the SystemState object and Settings objects are
    /// provided by storing them as shared pointers, which are initialized here.
    /// Memory on device is also allocated for several freeze out arrays used
    /// as cache.
    /// @param settingsPtr_in A shared pointer to the Settings object.
    /// @param systemPtr_in A shared pointer to the SystemState object.
    void initialize( std::shared_ptr<Settings> settingsPtr_in,
                     std::shared_ptr<SystemState<D>> systemPtr_in)
    {
      settingsPtr = settingsPtr_in;
      systemPtr = systemPtr_in;
      frz1 = freeze_array("frz1",     systemPtr->cabana_particles.size());
      frz2 = freeze_array("frz2",     systemPtr->cabana_particles.size());
      fback = freeze_array("fback",   systemPtr->cabana_particles.size());
      fback2 = freeze_array("fback2", systemPtr->cabana_particles.size());
      fback3 = freeze_array("fback3", systemPtr->cabana_particles.size());
      fback4 = freeze_array("fback4", systemPtr->cabana_particles.size());
      results = freeze_results("results", systemPtr->cabana_particles.size());
      mask = mask_type("mask",        systemPtr->cabana_particles.size());
    };

  public:
    freeze_results results; ///< freeze out results

    /// @brief Freeze-out procedure performed in the first time step
    /// @details This function is called in the first time step to set up the
    /// freeze-out snapshots at the first time step. This will fill the frz2
    /// array with the current state of the particles.
    void freezeout_first_time_step(){

      taupp = systemPtr->t;
      frzc  = 1;

      //Setup slices for frz2
      FRZ_VIEW(frz2_, frz2);
      CREATE_VIEW(device_, systemPtr->cabana_particles);
      auto setup_frz2_first_timestep = KOKKOS_CLASS_LAMBDA(const int is, const int ia){
        for (int idir=0;idir<D+1;++idir)
        for (int jdir=0;jdir<D+1;++jdir)
            frz2_shear.access(is, ia, idir, jdir) = device_hydro_spacetime_matrix.access(is, ia, ccake::hydro_info::shv, idir, jdir);
        for (int idir=0;idir<D;++idir){
          frz2_r.access(is, ia, idir) = device_position.access(is, ia, idir);
          frz2_u.access(is, ia, idir) = device_hydro_vector.access(is, ia, ccake::hydro_info::u, idir);
          frz2_gradP.access(is, ia, idir) = device_hydro_vector.access(is, ia, ccake::hydro_info::gradP, idir);
          frz2_gradE.access(is, ia, idir) = device_hydro_vector.access(is, ia, ccake::hydro_info::gradE, idir);
        }
        frz2_s.access(is, ia)    = device_thermo.access(is, ia, ccake::thermo_info::s);
        frz2_e.access(is, ia)    = device_thermo.access(is, ia, ccake::thermo_info::e);
        frz2_rhoB.access(is, ia) = device_thermo.access(is, ia, ccake::thermo_info::rhoB);
        frz2_rhoS.access(is, ia) = device_thermo.access(is, ia, ccake::thermo_info::rhoS);
        frz2_rhoQ.access(is, ia) = device_thermo.access(is, ia, ccake::thermo_info::rhoQ);
        frz2_T.access(is, ia)    = device_thermo.access(is, ia, ccake::thermo_info::T);
        frz2_muB.access(is, ia)  = device_thermo.access(is, ia, ccake::thermo_info::muB);
        frz2_muS.access(is, ia)  = device_thermo.access(is, ia, ccake::thermo_info::muS);
        frz2_muQ.access(is, ia)  = device_thermo.access(is, ia, ccake::thermo_info::muQ);
        frz2_theta.access(is, ia) = device_hydro_scalar.access(is, ia, ccake::hydro_info::div_u) + device_hydro_scalar.access(is, ia, ccake::hydro_info::gamma)/taupp;
        frz2_bulk.access(is, ia) = device_hydro_scalar.access(is, ia, ccake::hydro_info::bigPI);
        frz2_sigma.access(is, ia) = device_hydro_scalar.access(is, ia, ccake::hydro_info::sigma);
        frz2_shear33.access(is, ia) = device_hydro_scalar.access(is, ia, ccake::hydro_info::shv33);
        frz2_inside.access(is, ia) = device_hydro_scalar.access(is, ia, ccake::hydro_info::inside);
        frz2_wfz.access(is, ia) = device_thermo.access(is, ia, ccake::thermo_info::w);
        frz2_cs2fz.access(is, ia) = device_thermo.access(is, ia, ccake::thermo_info::cs2);
      };
      auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());
      Cabana::simd_parallel_for(simd_policy, setup_frz2_first_timestep, "setup_frz2_first_timestep");

    };

    /// @brief Freeze-out procedure performed in the second time step
    /// @details This function is called in the second time step to set up the
    /// freeze-out snapshots at the second time step. This will fill the frz1
    /// array with the current state of the particles.
    void freezeout_second_time_step(){
      taup = systemPtr->t;
      frzc = 2;
      //Setup slices for frz1
      FRZ_VIEW(frz1_, frz1);

      CREATE_VIEW(device_, systemPtr->cabana_particles);
      auto setup_frz1_first_timestep = KOKKOS_CLASS_LAMBDA(const int is, const int ia){
        for (int idir=0;idir<D+1;++idir)
        for (int jdir=0;jdir<D+1;++jdir)
            frz1_shear.access(is, ia, idir, jdir) = device_hydro_spacetime_matrix.access(is, ia, ccake::hydro_info::shv, idir, jdir);
        for (int idir=0;idir<D;++idir){
          frz1_r.access(is, ia, idir) = device_position.access(is, ia, idir);
          frz1_u.access(is, ia, idir) = device_hydro_vector.access(is, ia, ccake::hydro_info::u, idir);
          frz1_gradP.access(is, ia, idir) = device_hydro_vector.access(is, ia, ccake::hydro_info::gradP, idir);
          frz1_gradE.access(is, ia, idir) = device_hydro_vector.access(is, ia, ccake::hydro_info::gradE, idir);
        }
        frz1_s.access(is, ia)    = device_thermo.access(is, ia, ccake::thermo_info::s);
        frz1_e.access(is, ia)    = device_thermo.access(is, ia, ccake::thermo_info::e);
        frz1_rhoB.access(is, ia) = device_thermo.access(is, ia, ccake::thermo_info::rhoB);
        frz1_rhoS.access(is, ia) = device_thermo.access(is, ia, ccake::thermo_info::rhoS);
        frz1_rhoQ.access(is, ia) = device_thermo.access(is, ia, ccake::thermo_info::rhoQ);
        frz1_T.access(is, ia)    = device_thermo.access(is, ia, ccake::thermo_info::T);
        frz1_muB.access(is, ia)  = device_thermo.access(is, ia, ccake::thermo_info::muB);
        frz1_muS.access(is, ia)  = device_thermo.access(is, ia, ccake::thermo_info::muS);
        frz1_muQ.access(is, ia)  = device_thermo.access(is, ia, ccake::thermo_info::muQ);
        frz1_theta.access(is, ia) = device_hydro_scalar.access(is, ia, ccake::hydro_info::div_u) + device_hydro_scalar.access(is, ia, ccake::hydro_info::gamma)/taupp;
        frz1_bulk.access(is, ia) = device_hydro_scalar.access(is, ia, ccake::hydro_info::bigPI);
        frz1_sigma.access(is, ia) = device_hydro_scalar.access(is, ia, ccake::hydro_info::sigma);
        frz1_shear33.access(is, ia) = device_hydro_scalar.access(is, ia, ccake::hydro_info::shv33);
        frz1_inside.access(is, ia) = device_hydro_scalar.access(is, ia, ccake::hydro_info::inside);
        frz1_wfz.access(is, ia) = device_thermo.access(is, ia, ccake::thermo_info::w);
        frz1_cs2fz.access(is, ia) = device_thermo.access(is, ia, ccake::thermo_info::cs2);
      };
      auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());
      Cabana::simd_parallel_for(simd_policy, setup_frz1_first_timestep, "setup_frz1_first_timestep");
    }

    /// @brief Freeze-out procedure performed from the third time step onwards
    /// @details This function is used to update the snapshots of the freeze-out
    /// procedure from the third time step onwards. However, before the data
    /// movement, a mask is filled which will be used by the shift_snapshot 
    /// function to decide which particles should not be copied from one 
    /// snapshot to another.
    /// @see shift_snapshot
    void freeze_out_evo(){
      CREATE_VIEW(device_, systemPtr->cabana_particles);
      FRZ_VIEW(fback_, fback);

      // Create mask for case where freeze < 4 && btrack <= 3 && btrack > 0
      auto condition = Cabana::slice<0>(mask);
      auto tag_not_frozen_few_neighbors = KOKKOS_LAMBDA(const int is, const int ia){
        condition.access(is,ia) = (device_freeze.access(is,ia)< 4 &&
                                   device_btrack.access(is,ia) <= 3 &&
                                   device_btrack.access(is,ia) > 0);
      };
      auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());
      Cabana::simd_parallel_for(simd_policy, tag_not_frozen_few_neighbors, "tag_not_frozen_few_neighbors");
      Kokkos::fence();
      //And shifts the snapshots for this case
      shift_snapshot(fback4, fback2);
      shift_snapshot(fback3, fback);
      Kokkos::fence();
      shift_snapshot(fback2, frz2);
      shift_snapshot(fback, frz1);
      Kokkos::fence();

      // Now, we deal with the case where there is no neighbor (btrack == 0) and a non-zero gradient in the pressure (gradP == 0)
      auto tag_not_frozen_no_neighbors_finite_gradP =  KOKKOS_LAMBDA(const int is, const int ia){
        condition.access(is,ia) = (device_freeze.access(is,ia) < 4 &&
                                   device_btrack.access(is,ia) == 0 &&
                                   fback_gradP.access(is,ia, 0) != 0);
      };
      Cabana::simd_parallel_for(simd_policy, tag_not_frozen_no_neighbors_finite_gradP, "tag_not_frozen_no_neighbors_finite_gradP");
      Kokkos::fence();
      shift_snapshot(frz2, fback2);
      shift_snapshot(frz1, fback);
      Kokkos::fence();
      //Otherwise, we deal with the case where there is no neighbor (btrack == 0) and a zero gradient in the pressure (gradP == 0)
      auto tag_not_frozen_no_neighbors_zero_gradP =  KOKKOS_LAMBDA(const int is, const int ia){
        condition.access(is,ia) = (device_freeze.access(is,ia) < 4 &&
                                   device_btrack.access(is,ia) == 0 &&
                                   fback_gradP.access(is,ia, 0) == 0);
      };
      Cabana::simd_parallel_for(simd_policy, tag_not_frozen_no_neighbors_zero_gradP, "tag_not_frozen_no_neighbors_zero_gradP");
      Kokkos::fence();
      shift_snapshot(frz2, fback4);
      shift_snapshot(frz1, fback3);
      Kokkos::fence();

      // Either way, if the particle is not frozen out and has no neighbors, we freeze it out
      auto freeze_no_neighbors = KOKKOS_LAMBDA(const int is, const int ia){
        if (device_freeze.access(is,ia) < 4 && device_btrack.access(is,ia) == 0){
          device_freeze.access(is,ia) = 4;
          device_btrack.access(is,ia) = -1;
        }
      };
      Cabana::simd_parallel_for(simd_policy, freeze_no_neighbors, "freeze_no_neighbors");
      Kokkos::fence();
    }

    /// @brief Copy the contents of one freeze out snapshot to another
    /// @details This function copies the contents of one freeze out snapshot to
    /// another. Particles which are tagged in the mask are not copied though.
    /// @param frz_dest Destination freeze out snapshot
    /// @param frz_src Source freeze out snapshot
    void shift_snapshot(freeze_array&  frz_dest, freeze_array& frz_src)
    {
      FRZ_VIEW(frz_dest_, frz_dest);
      FRZ_VIEW(frz_src_, frz_src);
      auto perform_shift = Cabana::slice<0>(mask);

      auto shift = KOKKOS_CLASS_LAMBDA(const int is, const int ia){
        if (perform_shift.access(is,ia)){
          for (int idir=0;idir<D+1;++idir)
          for (int jdir=0;jdir<D+1;++jdir)
            frz_dest_shear.access(is, ia, idir, jdir) = frz_src_shear.access(is, ia, idir, jdir);
          for (int idir=0;idir<D;++idir){
            frz_dest_r.access(is, ia, idir) = frz_src_r.access(is, ia, idir);
            frz_dest_u.access(is, ia, idir) = frz_src_u.access(is, ia, idir);
            frz_dest_gradP.access(is, ia, idir) = frz_src_gradP.access(is, ia, idir);
            frz_dest_gradE.access(is, ia, idir) = frz_src_gradE.access(is, ia, idir);
          }
          frz_dest_t.access(is, ia) = frz_src_t.access(is, ia);
          frz_dest_s.access(is, ia) = frz_src_s.access(is, ia);
          frz_dest_e.access(is, ia) = frz_src_e.access(is, ia);
          frz_dest_rhoB.access(is, ia) = frz_src_rhoB.access(is, ia);
          frz_dest_rhoS.access(is, ia) = frz_src_rhoS.access(is, ia);
          frz_dest_rhoQ.access(is, ia) = frz_src_rhoQ.access(is, ia);
          frz_dest_T.access(is, ia) = frz_src_T.access(is, ia);
          frz_dest_muB.access(is, ia) = frz_src_muB.access(is, ia);
          frz_dest_muS.access(is, ia) = frz_src_muS.access(is, ia);
          frz_dest_muQ.access(is, ia) = frz_src_muQ.access(is, ia);
          frz_dest_theta.access(is, ia) = frz_src_theta.access(is, ia);
          frz_dest_bulk.access(is, ia) = frz_src_bulk.access(is, ia);
          frz_dest_sigma.access(is, ia) = frz_src_sigma.access(is, ia);
          frz_dest_shear33.access(is, ia) = frz_src_shear33.access(is, ia);
          frz_dest_inside.access(is, ia) = frz_src_inside.access(is, ia);
          frz_dest_wfz.access(is, ia) = frz_src_wfz.access(is, ia);
          frz_dest_cs2fz.access(is, ia) = frz_src_cs2fz.access(is, ia);
        }
      };
      auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());
      Cabana::simd_parallel_for(simd_policy, shift, "shift");
    }

    //==========================================================================
    //Constructor
    FreezeOut() = delete;
    /// @brief Default constructor
    /// @see initialize
    /// @param settingsPtr_in 
    /// @param systemPtr_in 
    FreezeOut( std::shared_ptr<Settings> settingsPtr_in,
               std::shared_ptr<SystemState<D>> systemPtr_in)
    {
      initialize( settingsPtr_in, systemPtr_in );
    };

    /// @brief Checks the freeze-out status of particles and updates the count of particles freezing out.
    /// @details The freeze-out procedure is as follows:
    /// - If the energy of a particle is below the freeze-out energy and
    /// the freeze-out has not yet begun (status 1), the particle is marked as 
    /// starting to freeze out by changing its status to 1.
    /// - If the particle started freezing out in the previous time step 
    /// (status 1), the following may happen
    ///   * If the particle has no neighbors or if its energy still below the 
    ///     freeze-out energy, the particle is marked as frozen out by updating
    ///     its status to 3.
    ///   * If the particle has neighbors and its energy has increased since the
    ///     previous time step, the particle is still a candidate to freeze out
    ///     and its status is kept as 1.
    ///   * If all the above is false, then the particle is not ready to freeze
    ///     out and its status is set to 0.
    /// @param count [in, out] The count of particles freezing out.
    void check_freeze_out_status( int& count )
    {
      double time_in = systemPtr->t;
      double efcheck = systemPtr->efcheck;
      auto frz2_t = Cabana::slice<FRZ_enum::t>(frz2);
      auto frz1_t = Cabana::slice<FRZ_enum::t>(frz1);
      auto frz1_e = Cabana::slice<FRZ_enum::e>(frz1);
      CREATE_VIEW(device_, systemPtr->cabana_particles);

      auto set_status = KOKKOS_LAMBDA(const int is,  const int ia){
        int freeze =  device_freeze.access(is,ia);
        double efchechk = device_efcheck.access(is,ia);
        double e = device_thermo.access(is, ia, ccake::thermo_info::e);

        switch (freeze)
        {
        case 0:                                   // if particle has not yet started freezing out
          if ( e <= efcheck ){                    // and if it's energy is below freeze out energy
            device_freeze.access(is,ia) = 1;         // then we start freeze-out
            frz2_t.access(is,ia) = time_in;
          }
          break;
        case 1:                                   // otherwise, if freeze-out has begun
          if ( device_btrack.access(is,ia) == -1 ){  // and if particle has no neighbors
            //count += 1;                           // then increment currently freezing out count
            device_freeze.access(is,ia) = 3;         // and freezes the particle.
            frz1_t.access(is,ia) = time_in;
          } else if ( e > frz1_e.access(is,ia) ){    // otherwise, if the particle's energy has increased has increased since the previous timestep
            device_freeze.access(is,ia) = 1;         // keep attempting freeze-out
            frz2_t.access(is,ia) = time_in;
          } else if( e <= efcheck ){              // otherwise, if e is below eFO
            //count += 1;                           // increment currently freezing out count
            device_freeze.access(is,ia) = 3;         // and freezes the particle.
            frz1_t.access(is,ia) = time_in;
          } else {                                // otherwise, particle is not ready to freeze-out
            device_freeze.access(is,ia)=0;
          }
          break;
        default:
          //Any other case, particle is already frozen out
          //Do nothing
          break;
        }
      };
      auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());
      Cabana::simd_parallel_for(simd_policy, set_status, "set_status");
      Kokkos::fence();
      //Count how many particles are freezing out
      count=0;
      auto count_freeze = KOKKOS_LAMBDA(const int& iparticle, int& count){
        count += (device_freeze(iparticle) == 3) ? 1 : 0;
      };
      Kokkos::parallel_reduce("count_freeze",
                              systemPtr->cabana_particles.size(), count_freeze, count);
      Kokkos::fence();
      systemPtr->number_part_fo += count;
      //Copy aux to host to make freeze report
      #ifdef DEBUG_SLOW
      auto h_particles = Cabana::create_mirror_view_and_copy(Kokkos::HostSpace(), systemPtr->cabana_particles);
      CREATE_VIEW(h_, h_particles);

      //LOOP OVER PARTICLES AND COUNT HOW MANY ARE FREEZING OUT
      int count_status0 = 0;
      int count_status1 = 0;
      int count_status3 = 0;
      int count_status4 = 0;
      for(int ipart=0; ipart<systemPtr->cabana_particles.size(); ++ipart ){
        count_status0 += (h_freeze(ipart) < 0.5) ? 1 : 0;
        count_status1 += (h_freeze(ipart) > 0.5 && h_freeze(ipart) < 1.5) ? 1 : 0;
        count_status3 += (h_freeze(ipart) > 2.5 && h_freeze(ipart) < 3.5) ? 1 : 0;
        count_status4 += (h_freeze(ipart) > 3.5) ? 1 : 0;
      }
      formatted_output::detail("Freeze energy...................: "+std::to_string(efcheck*constants::hbarc_GeVfm));
      formatted_output::detail("Number of particles freezing out: "+std::to_string(count));
      formatted_output::detail("Number of particles status 0....: "+std::to_string(count_status0));
      formatted_output::detail("Number of particles status 1....: "+std::to_string(count_status1));
      formatted_output::detail("Number of particles status 3....: "+std::to_string(count_status3));
      formatted_output::detail("Number of particles status 4....: "+std::to_string(count_status4));

      #endif
      //Last, we clear the tag for printing
      FRZ_RESULTS_VIEW(results_, results)
      auto clear_print = KOKKOS_LAMBDA(const int& iparticle){
        results_print(iparticle) = false;
      };
      Kokkos::parallel_for("clear_print", systemPtr->cabana_particles.size(), clear_print);
      Kokkos::fence();
    }

    /// @brief Shell function to control the freeze-out procedure.
    /// @details This function controls the freeze-out procedure. It controls
    /// the caches of freeze-out snapshots, and calls the interpolator function
    /// to interpolate the freeze-out quantities between the snapshots.
    /// @param curfrz
    void bsqsvfreezeout(int curfrz)
    {
      if (frzc==0) // true in first timestep only
      {
        #ifdef DEBUG
        formatted_output::detail("First timestep freeze-out");
        #endif
        freezeout_first_time_step();
      }
      else if (frzc==1) // true in second timestep only
      {
        #ifdef DEBUG
        formatted_output::detail("Second timestep freeze-out");
        #endif
        freezeout_second_time_step();
        if ( curfrz > 0 )
          bsqsvinterpolate( curfrz );
        else
          cf = 0;
      }
      else // true in third timestep and afterward
      {
        freeze_out_evo();
        tau = systemPtr->t;

        if ( curfrz > 0 )
          bsqsvinterpolate( curfrz );
        else
          cf = 0;

        //sets up the variables for the next time step
        FRZ_VIEW(frz1_, frz1);
        CREATE_VIEW(device_, systemPtr->cabana_particles);
        //Setup next step_step
        auto condition = Cabana::slice<0>(mask);
        auto tag_all_to_swap = KOKKOS_LAMBDA(const int is, const int ia){
          condition.access(is,ia) = device_freeze.access(is,ia) == 4;
        };
        auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());
        Cabana::simd_parallel_for(simd_policy, tag_all_to_swap, "tag_all_to_swap");
        Kokkos::fence();
        shift_snapshot(frz2, frz1);
        Kokkos::fence();

        auto setup_next_step = KOKKOS_CLASS_LAMBDA(const int is, const int ia){
          for (int idir=0;idir<D+1;++idir)
          for (int jdir=0;jdir<D+1;++jdir)
            frz1_shear.access(is, ia, idir, jdir) = device_hydro_spacetime_matrix.access(is, ia, ccake::hydro_info::shv, idir, jdir);
          for (int idir=0;idir<D;++idir){
            frz1_r.access(is, ia, idir) = device_position.access(is, ia, idir);
            frz1_u.access(is, ia, idir) = device_hydro_vector.access(is, ia, ccake::hydro_info::u, idir);
            frz1_gradP.access(is, ia, idir) = device_hydro_vector.access(is, ia, ccake::hydro_info::gradP, idir);
            frz1_gradE.access(is, ia, idir) = device_hydro_vector.access(is, ia, ccake::hydro_info::gradE, idir);

          }
          frz1_s.access(is, ia) = device_thermo.access(is, ia, ccake::thermo_info::s);
          frz1_e.access(is, ia) = device_thermo.access(is, ia, ccake::thermo_info::e);
          frz1_rhoB.access(is, ia) = device_thermo.access(is, ia, ccake::thermo_info::rhoB);
          frz1_rhoS.access(is, ia) = device_thermo.access(is, ia, ccake::thermo_info::rhoS);
          frz1_rhoQ.access(is, ia) = device_thermo.access(is, ia, ccake::thermo_info::rhoQ);
          frz1_T.access(is, ia) = device_thermo.access(is, ia, ccake::thermo_info::T);
          frz1_muB.access(is, ia) = device_thermo.access(is, ia, ccake::thermo_info::muB);
          frz1_muS.access(is, ia) = device_thermo.access(is, ia, ccake::thermo_info::muS);
          frz1_muQ.access(is, ia) = device_thermo.access(is, ia, ccake::thermo_info::muQ);
          frz1_theta.access(is, ia) = device_hydro_scalar.access(is, ia, ccake::hydro_info::div_u) + device_hydro_scalar.access(is, ia, ccake::hydro_info::gamma)/tau;
          frz1_bulk.access(is, ia) = device_hydro_scalar.access(is, ia, ccake::hydro_info::bigPI);
          frz1_sigma.access(is, ia) = device_hydro_scalar.access(is, ia, ccake::hydro_info::sigma);
          frz1_shear33.access(is, ia) = device_hydro_scalar.access(is, ia, ccake::hydro_info::shv33);
          frz1_inside.access(is, ia) = device_hydro_scalar.access(is, ia, ccake::hydro_info::inside);
          frz1_wfz.access(is, ia) = device_thermo.access(is, ia, ccake::thermo_info::w);
          frz1_cs2fz.access(is, ia) = device_thermo.access(is, ia, ccake::thermo_info::cs2);
        };
        Cabana::simd_parallel_for(simd_policy, setup_next_step, "setup_next_step");

        taupp = taup;
        taup  = tau;
      }
    }

    /// @brief Interpolates the freeze-out quantities between the snapshots
    /// @details This function interpolates the freeze-out quantities between the
    /// snapshots. It decides which snapshot to use for the interpolation based
    /// on the distance of the particle to the freeze-out energy. The function
    /// then interpolates the freeze-out quantities between the snapshots. Only
    /// particles which are currently freezing out (status = 4) are interpolated.
    /// @param curfrz
    void bsqsvinterpolate(int curfrz)
    {

      double efcheck = systemPtr->efcheck;

      FRZ_VIEW(frz1_, frz1);
      FRZ_VIEW(frz2_, frz2);
      FRZ_RESULTS_VIEW(results_, results);

      #ifdef DEBUG_SLOW
      auto h_particles = Cabana::create_mirror_view_and_copy(Kokkos::HostSpace(), systemPtr->cabana_particles);
      CREATE_VIEW(h_, h_particles);

      //LOOP OVER PARTICLES AND COUNT HOW MANY ARE FREEZING OUT
      int count_status0 = 0;
      int count_status1 = 0;
      int count_status3 = 0;
      int count_status4 = 0;
      for(int ipart=0; ipart<systemPtr->cabana_particles.size(); ++ipart ){
        count_status0 += (h_freeze(ipart) < 0.5) ? 1 : 0;
        count_status1 += (h_freeze(ipart) > 0.5 && h_freeze(ipart) < 1.5) ? 1 : 0;
        count_status3 += (h_freeze(ipart) > 2.5 && h_freeze(ipart) < 3.5) ? 1 : 0;
        count_status4 += (h_freeze(ipart) > 3.5) ? 1 : 0;
      }
      formatted_output::detail("Number of particles status 3 before interpolate....: "+std::to_string(count_status3));
      formatted_output::detail("Number of particles status 4 before interpolate....: "+std::to_string(count_status4));
      #endif
      //Last, we clear the tag for printing
      auto clear_print = KOKKOS_LAMBDA(const int& iparticle){
        results_print(iparticle) = false;
      };

      double dt = settingsPtr->dt;
      // for all CURRENTLY FREEZING OUT PARTICLES
      CREATE_VIEW(device_, systemPtr->cabana_particles);
      auto interpolate = KOKKOS_CLASS_LAMBDA(const int is, const int ia){
        if (device_freeze.access(is,ia) == 3){ //If particles are tagged to freeze out

          // decide whether the particle was closer to freeze out at the previous
          // timestep or the one before that
          int swit = ( abs( frz1_e.access(is, ia) - efcheck ) < abs( frz2_e.access(is, ia) - efcheck ) ) ? 1 : 2;

          double sigsub = 0.0, inside = 0.0;
          milne::Vector<double,D> gradPsub, gradEsub;
          if ( swit == 1 )  // if particle was closer to freeze-out at last timestep
          {
            // if particle had neighbors, use previous timestep otherwise go back two timesteps
            results_tlist.access(is,ia) = device_btrack.access(is,ia) != -1 ? taup : taup - dt;
            for (int idir=0; idir<D+1; ++idir)
            for (int jdir=0; jdir<D+1; ++jdir)
              results_shearsub.access(is,ia,idir,jdir) = frz1_shear.access(is,ia,idir,jdir);
            for (int idir=0; idir<D; ++idir){
              results_rsub.access(is,ia,idir) = frz1_r.access(is,ia,idir);
              results_uout.access(is,ia,idir) = frz1_u.access(is,ia,idir);
              gradPsub(idir) = frz1_gradP.access(is,ia,idir);
              gradEsub(idir) = frz1_gradE.access(is,ia,idir);
            }
            results_bulksub.access(is,ia) = frz1_bulk.access(is,ia);
            results_thetasub.access(is,ia) = frz1_theta.access(is,ia);
            inside = frz1_inside.access(is,ia);
            sigsub = frz1_sigma.access(is,ia);
            results_shear33sub.access(is,ia) = frz1_shear33.access(is,ia);
            results_Efluc.access(is,ia) = frz1_e.access(is,ia);
            results_Tfluc.access(is,ia) = frz1_T.access(is,ia);
            results_muBfluc.access(is,ia) = frz1_muB.access(is,ia);
            results_muSfluc.access(is,ia) = frz1_muS.access(is,ia);
            results_muQfluc.access(is,ia) = frz1_muQ.access(is,ia);
            results_sFO.access(is,ia) = frz1_s.access(is,ia);
            results_cs2fzfluc.access(is,ia) = frz1_cs2fz.access(is,ia);
            results_wfzfluc.access(is,ia) = frz1_wfz.access(is,ia);
          } else if (swit == 2){
            // if particle had neighbors, use previous timestep otherwise go back two timesteps
            results_tlist.access(is,ia) = device_btrack.access(is,ia) != -1 ? taupp : taupp - dt;
            for (int idir=0; idir<D+1; ++idir)
            for (int jdir=0; jdir<D+1; ++jdir)
              results_shearsub.access(is,ia,idir,jdir) = frz2_shear.access(is,ia,idir,jdir);
            for (int idir=0; idir<D; ++idir){
              results_rsub.access(is,ia,idir) = frz2_r.access(is,ia,idir);
              results_uout.access(is,ia,idir) = frz2_u.access(is,ia,idir);
              gradPsub(idir) = frz2_gradP.access(is,ia,idir);
              gradEsub(idir) = frz2_gradE.access(is,ia,idir);
            }
            results_bulksub.access(is,ia) = frz2_bulk.access(is,ia);
            results_thetasub.access(is,ia) = frz2_theta.access(is,ia);
            inside = frz2_inside.access(is,ia);
            sigsub = frz2_sigma.access(is,ia);
            results_shear33sub.access(is,ia) = frz2_shear33.access(is,ia);
            results_Efluc.access(is,ia) = frz2_e.access(is,ia);
            results_Tfluc.access(is,ia) = frz2_T.access(is,ia);
            results_muBfluc.access(is,ia) = frz2_muB.access(is,ia);
            results_muSfluc.access(is,ia) = frz2_muS.access(is,ia);
            results_muQfluc.access(is,ia) = frz2_muQ.access(is,ia);
            results_sFO.access(is,ia) = frz2_s.access(is,ia);
            results_cs2fzfluc.access(is,ia) = frz2_cs2fz.access(is,ia);
            results_wfzfluc.access(is,ia) = frz2_wfz.access(is,ia);
          }

          // COMPUTE NORMALS AFTER THIS POINT
          double norm2 = 0.0;
          for (int idir=0; idir<D; ++idir) norm2 += results_uout.access(is,ia,idir)*results_uout.access(is,ia,idir);
          results_gsub.access(is,ia) = Kokkos::sqrt( norm2 + 1 );
          sigsub /= results_gsub.access(is,ia)*results_tlist.access(is,ia);
          results_swsub.access(is,ia) = device_norm_spec.access(is, ia, ccake::densities_info::s)/sigsub;
          for(int idir=0;idir<D;++idir){
            results_divT.access(is,ia,idir) = (1.0/results_sFO.access(is,ia))*gradPsub(idir);//\partial_\mu T= \partial_\mu P/s.
            results_divP.access(is, ia, idir) = gradPsub(idir);
            results_divE.access(is, ia, idir) = gradEsub(idir);
          }
          double cs2 = results_cs2fzfluc.access(is,ia);
          double w = results_wfzfluc.access(is,ia);
          double bulk = results_bulksub.access(is,ia);
          double theta = results_thetasub.access(is,ia);

          double inner = 0;
          for (int idir=0; idir<D; ++idir) inner += results_uout.access(is,ia,idir)*gradPsub(idir);
          results_divPpress.access(is,ia) = - (1.0/results_gsub.access(is,ia))*( cs2*(w+bulk)*theta - cs2*inside + inner ); //\partial_\tau P: Eq C3 in arXiv.1406.3333
          results_divTtemp.access(is,ia) = results_divPpress.access(is,ia)/results_sFO.access(is,ia);
          inner = 0;
          for (int idir=0; idir<D; ++idir) inner += results_uout.access(is,ia,idir)*gradEsub(idir);
          results_divEener.access(is,ia) = - (1.0/results_gsub.access(is,ia))*( (w+bulk)*theta - inside + inner); //De -> \partial_\tau e: Eq C2 in arXiv.1406.3333
          //THIS NEEDS TO BE RESET
          norm2 = 0;
          for (int idir=0; idir<D; ++idir) norm2 += results_divT.access(is,ia,idir)*results_divT.access(is,ia,idir);
          double insub = results_divTtemp.access(is,ia)*results_divTtemp.access(is,ia) - norm2;
          double norm  = -Kokkos::sqrt(Kokkos::fabs(insub));
          norm = Kokkos::fabs(norm) > 1.E-14 ? norm : 1.E-14;
          results_divTtemp.access(is,ia) /= norm;
          for (int idir=0; idir<D; ++idir) results_divT.access(is,ia,idir) /= norm;

          norm2 = 0;
          for (int idir=0; idir<D; ++idir) norm2 += results_divE.access(is,ia,idir)*results_divE.access(is,ia,idir);
          double insubE = results_divEener.access(is,ia)*results_divEener.access(is,ia) - norm2;
          double normE  = -Kokkos::sqrt(Kokkos::fabs(insubE));
          normE = Kokkos::fabs(normE) > 1.E-14 ? normE : 1.E-14;
          results_divEener.access(is,ia) /= normE;
          for(int idir=0;idir<D;++idir) results_divE.access(is,ia,idir) /= normE;

          results_sFO.access(is,ia) *= Kokkos::pow(results_Tfluc.access(is,ia)*constants::hbarc_GeVfm, 3);
          results_Tfluc.access(is,ia) *= constants::hbarc_GeVfm;

          // Lastly, we tag the particle as frozen and to be printed
          device_freeze.access(is,ia) = 4;
          results_print.access(is,ia) = true;
        } //end if particle is frozen
      }; //end interpolate lambda
      auto simd_policy = Cabana::SimdPolicy<VECTOR_LENGTH,ExecutionSpace>(0, systemPtr->cabana_particles.size());
      Cabana::simd_parallel_for(simd_policy, interpolate, "interpolate");
      Kokkos::fence();
      #ifdef DEBUG_SLOW
      Cabana::deep_copy(h_freeze, device_freeze);
      count_status3 = 0;
      count_status4 = 0;
      for(int ipart=0; ipart<systemPtr->cabana_particles.size(); ++ipart ){
        count_status3 += (h_freeze(ipart) > 2.5 && h_freeze(ipart) < 3.5) ? 1 : 0;
        count_status4 += (h_freeze(ipart) > 3.5) ? 1 : 0;
      }
      formatted_output::detail("Number of particles status 3 after interpolate....: "+std::to_string(count_status3));
      formatted_output::detail("Number of particles status 4 after interpolate....: "+std::to_string(count_status4));
      #endif
      cf = curfrz;
    }
};
}

#endif