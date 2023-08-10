#ifndef THERMODYNAMIC_INFO_H
#define THERMODYNAMIC_INFO_H

namespace ccake{

/// \brief struct to hold thermodynamic info of the particle.

/// NOTE: Beware these information needs to be passed to the cabana AoSoA. To this
/// end, we define helper macros at the bottom of this file. If one wants to
/// add a new member to this struct, one needs to add it to the macros as well.
/// Also, remmember to pass new datamembers to cabana AoSoA in
/// `SystemState::allocate_cababa_particles()`
struct thermodynamic_info
{
  // (T,mu_i) coordinates depend on which EoS was used!
  string eos_name = "";

  double T    = 0.0, muB  = 0.0, muS  = 0.0, muQ  = 0.0;
  double e    = 0.0, s    = 0.0, rhoB = 0.0, rhoS = 0.0, rhoQ = 0.0,
         p    = 0.0, cs2  = 0.0, w    = 0.0, A    = 0.0;
  double dwds = 0.0, dwdB = 0.0, dwdS = 0.0, dwdQ = 0.0;

};

namespace thermo_info{
  enum thermo_scalar_info{
    T,
    muB,
    muS,
    muQ,
    e,
    s,
    rhoB,
    rhoS,
    rhoQ,
    p,
    cs2,
    w,
    A,
    dwds,
    dwdB,
    dwdS,
    dwdQ,
    NUM_THERMO_INFO
  };
#define THERMO_SCALAR_INFO double[ccake::thermo_info::NUM_THERMO_INFO]
}}

#endif