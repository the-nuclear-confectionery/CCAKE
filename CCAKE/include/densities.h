#ifndef DENSITIES_H
#define DENSITIES_H


/// \brief struct to hold densities
/// \details This struct is used to hold the densities carried by an SPH particle.

/// NOTE: Beware these information needs to be passed to the cabana AoSoA. To this
/// end, we define helper macros at the bottom of this file. If one wants to
/// add a new member to this struct, one needs to add it to the macros as well.
/// Also, remmember to pass new datamembers to cabana AoSoA in
/// `SystemState::allocate_cababa_particles()`
namespace ccake
{
struct densities
{
  double e    = 0.0;
  double s    = 0.0;
  double rhoB = 0.0;
  double rhoS = 0.0;
  double rhoQ = 0.0;
};

namespace densities_info{
  enum density_info{
    e,
    s,
    rhoB,
    rhoS,
    rhoQ,
    NUM_DENSITIES
  };
  #define DENSITY_INFO double[ccake::densities_info::NUM_DENSITIES]
}}

#endif