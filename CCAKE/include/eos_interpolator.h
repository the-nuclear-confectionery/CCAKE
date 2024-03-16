#ifndef EOS_INTERPOLATOR_H
#define EOS_INTERPOLATOR_H

#include <filesystem>
#include <iostream>
#include <cmath>

#include <Kokkos_Core.hpp>
#include <Cabana_Core.hpp>

#include "H5Cpp.h"

#include "utilities.h"
#include "particle.h"


namespace fs = std::filesystem;

/// @brief Class for interpolating the EoS table
/// Assumes that the EoS table is in the format:
/// The first 4 lines are headers, containing the limits of the grid in each dimension
/// There is a 5th header line with the names of the variables
/// The rest of the lines are space separated values of the variables
/// The table has as independent variables: s, rhoB, rhoQ, rhoS. These are laid
/// out in a regular grid. The slowest changing variable is rhoS, followed by rhoQ,
/// followed by rhoB, followed by s.

namespace ccake {

using eos_thermo = Kokkos::View<const double****, DeviceType, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;
using eos_thermo_nonconst = Kokkos::View<double****, DeviceType>;
using eos_thermo_kind = Kokkos::View<const int****, DeviceType, Kokkos::MemoryTraits<Kokkos::RandomAccess>>;

namespace eos_variables{
enum eos_variables_enum{
  T, muB, muQ, muS,
  e, p, cs2,
  dw_ds, dw_dB, dw_dQ, dw_dS,
  NUM_EOS_VARIABLES
  };
}

class EoS_Interpolator{

private:
  //eos_thermo *eos_vars[eos_variables::NUM_EOS_VARIABLES];
  eos_thermo_kind *eos_type;
  eos_thermo_nonconst eos_vars[eos_variables::NUM_EOS_VARIABLES];
  //eos_thermo_nonconst T, muB, muQ, muS, e, p, cs2, dw_ds, dw_dB, dw_dQ, dw_dS;

  //Parameters of the grid in each dimension
  double s_min, s_max, ds,
         B_min, B_max, dB,
         S_min, S_max, dS,
         Q_min, Q_max, dQ;
  int Ns, NB, NS, NQ;

  KOKKOS_FUNCTION
  double interpolate4D(int idx[], double pos[], int ivar) const;

public:
  EoS_Interpolator() = delete;
  EoS_Interpolator(fs::path path_to_eos_table);
  ~EoS_Interpolator() = default;

  void fill_thermodynamics(Cabana::AoSoA<CabanaParticle, DeviceType, VECTOR_LENGTH> &particles, const double t);

};

} // namespace ccake
#endif