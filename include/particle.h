#ifndef PARTICLE_H
#define PARTICLE_H

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <Cabana_Core.hpp>

#include "constants.h"
#include "densities.h"
#include "hydrodynamic_info.h"
#include "matrix.h"
#include "settings.h"
#include "thermodynamic_info.h"
#include "vector.h"

namespace ccake {
///=============================================================================
///\brief Particle definitions for use in Cabana library
///\details This is an alias for Cabana::AoSoA, which is the data structure.
/// This allows for efficient memory access and data management in the library.

///IMPORTANT: For optimal performance, the particle data should be listed from
///larger size to smaller size (i.e.,double[D], double, float, int, etc.). See
///cabana documentation on AoSoA for more details. Notice also that usagge of
///non arithmetic types (e.g., Vector<double,D>) requires the use of a Cabana
///member type (see below).
///
///NOTE: Only data types that are arithmetic types (i.e., int, float, double,
///long, etc.) can be used as template parameters for Cabana::AoSoA. Arrays of
///arithmetic types are allowed as well. However, classes and structs are
///disallowed.
///
/// NOTE: Beware that some datamembers needs to be passed to the cabana AoSoA.
/// These should be specified in the `CabanaParticle` alias defined below. One
/// should also add it to the enum `particle_data` defined below. Order matters.
/// Please, also create the necessary slice in the CREATE_VIEW macro.
/// Also, remmember to pass new datamembers to cabana AoSoA in
/// `SystemState::allocate_cababa_particles()`
///
///TODO: I have almost certain that we do not need to load all of these data
/// members into the device. We should probably only load the data members that
/// are needed for the evolution. A cleanup may be necessary in the future.

using CabanaParticle = Cabana::MemberTypes<
                                        HYDRO_SPACE_MATRIX_INFO, //72 doubles
                                        HYDRO_VECTOR_INFO,       //51 doubles
                                        HYDRO_SHEAR_AUX_VECTOR_INFO, // 32 doubles
                                        HYDRO_SCALAR_INFO,       //31 doubles
                                        HYDRO_SHEAR_AUX_MATRIX_INFO, // 21 doubles
                                        THERMO_SCALAR_INFO,      //17 doubles
                                        HYDRO_SPACETIME_MATRIX_INFO, // 16 doubles
                                        DENSITY_INFO,            // 5 doubles - INPUT
                                        DENSITY_INFO,            // 5 doubles - SMOOTHED
                                        DENSITY_INFO,            // 5 doubles - SPECIFIC
                                        DENSITY_INFO,            // 5 doubles - D_DT_SPEC
                                        DENSITY_INFO,            // 5 doubles - NORM_SPEC
                                        double[3], // position
                                        double,     // efcheck
                                        double,     // contribution_to_total_E
                                        double,     // contribution_to_total_dEz
                                        double,     // contribution_to_total_Ez
                                        int,        // particle id
                                        int,        // btrack
                                        int         // freeze
                                        >;
namespace particle_info
{
enum particle_data{
  hydro_space_matrix_info,
  hydro_vector_info,
  hydro_shear_aux_vector_info,
  hydro_scalar_info,
  hydro_shear_aux_matrix_info,
  thermo_scalar_info,
  hydro_spacetime_matrix_info,
  input_density,
  smoothed_density,
  specific_density,
  d_dt_spec_density,
  norm_spec_density,
  position,
  efcheck,
  contribution_to_total_E,
  contribution_to_total_dEz,
  contribution_to_total_Ez,
  id,
  btrack,
  freeze
};
}


//Macro to concatenate variable names
#ifndef CONCAT
  #define CONCAT_(a, b) a##b
  #define CONCAT(a, b) CONCAT_(a, b)
#endif

//Helper macros for creating views
#define CREATE_VIEW(prefix,cabana_aosoa) \
  auto CONCAT(prefix,hydro_space_matrix) = Cabana::slice<ccake::particle_info::hydro_space_matrix_info>(cabana_aosoa); \
  auto CONCAT(prefix,hydro_vector) = Cabana::slice<ccake::particle_info::hydro_vector_info>(cabana_aosoa); \
  auto CONCAT(prefix,hydro_shear_aux_matrix) = Cabana::slice<ccake::particle_info::hydro_shear_aux_matrix_info>(cabana_aosoa); \
  auto CONCAT(prefix,hydro_scalar) = Cabana::slice<ccake::particle_info::hydro_scalar_info>(cabana_aosoa); \
  auto CONCAT(prefix,hydro_shear_aux_vector) = Cabana::slice<ccake::particle_info::hydro_shear_aux_vector_info>(cabana_aosoa); \
  auto CONCAT(prefix,thermo) = Cabana::slice<ccake::particle_info::thermo_scalar_info>(cabana_aosoa); \
  auto CONCAT(prefix,hydro_spacetime_matrix) = Cabana::slice<ccake::particle_info::hydro_spacetime_matrix_info>(cabana_aosoa); \
  auto CONCAT(prefix,input) = Cabana::slice<ccake::particle_info::input_density>(cabana_aosoa); \
  auto CONCAT(prefix,smoothed) = Cabana::slice<ccake::particle_info::smoothed_density>(cabana_aosoa); \
  auto CONCAT(prefix,specific_density) = Cabana::slice<ccake::particle_info::specific_density>(cabana_aosoa); \
  auto CONCAT(prefix,d_dt_spec) = Cabana::slice<ccake::particle_info::d_dt_spec_density>(cabana_aosoa); \
  auto CONCAT(prefix,norm_spec) = Cabana::slice<ccake::particle_info::norm_spec_density>(cabana_aosoa); \
  auto CONCAT(prefix,position) = Cabana::slice<ccake::particle_info::position>(cabana_aosoa); \
  auto CONCAT(prefix,efcheck) = Cabana::slice<ccake::particle_info::efcheck>(cabana_aosoa); \
  auto CONCAT(prefix,contribution_to_total_E) = Cabana::slice<ccake::particle_info::contribution_to_total_E>(cabana_aosoa); \
  auto CONCAT(prefix,contribution_to_total_dEz) = Cabana::slice<ccake::particle_info::contribution_to_total_dEz>(cabana_aosoa); \
  auto CONCAT(prefix,contribution_to_total_Ez) = Cabana::slice<ccake::particle_info::contribution_to_total_Ez>(cabana_aosoa); \
  auto CONCAT(prefix,id) = Cabana::slice<ccake::particle_info::id>(cabana_aosoa); \
  auto CONCAT(prefix,btrack) = Cabana::slice<ccake::particle_info::btrack>(cabana_aosoa); \
  auto CONCAT(prefix,freeze) = Cabana::slice<ccake::particle_info::freeze>(cabana_aosoa);

template<unsigned int D>
class Particle
{
  private:

    static constexpr int VERBOSE = 0;

    Settings * settingsPtr   = nullptr;

  public:

    //==========================================================================
    // MEMBERS
    bool print_this_particle         = false;

    int ID                           = -1; ///< Particle ID
    int btrack                       = 0;  ///< Backtrack flag (number of neighbors)
    int Freeze                       = 0;  ///< Freeze flag

    double efcheck                   = 0.0;
    double contribution_to_total_E   = 0.0; ///< Contribution to total energy
    double contribution_to_total_dEz = 0.0;
    double contribution_to_total_Ez  = 0.0;

    Vector<double,D> r;       ///< Particle position

    //==========================================================================
    // different combinations of densities
    densities input     = {}; ///< These densities are read in initially and have physical units ~1/fm^3
    densities smoothed  = {}; ///< These are the smoothed (propagated) densities which have units ~/1/fm^2;
    densities specific  = {}; ///< These are the densities "per particle" which are effectively dimensionless;
    densities d_dt_spec = {}; ///< These are the time derivatives of the specific densities
    densities norm_spec = {}; ///< Gives the normalizations of the specific densities. Can choose different values for different densities by convenience;

    hydrodynamic_info<D> hydro  = {}; ///< Hydrodynamic information
    thermodynamic_info thermo = {}; ///< Thermodynamic information

    // Constructors and destructors.
    Particle();

    // copy-constructor
    Particle( const Particle& p );
   ~Particle(){}

    bool operator==( const Particle & ) const;

    // use this to set equation of state object before creating particles
    void set_SettingsPtr(Settings * settingsPtr_in) { settingsPtr = settingsPtr_in; }

    void evaluate_time_derivatives( double t );

    //==========================================================================
    // getter functions for thermodynamic information
    double T()    { return thermo.T;    }
    double muB()  { return thermo.muB;  }
    double muS()  { return thermo.muS;  }
    double muQ()  { return thermo.muQ;  }

    double p()    { return thermo.p;    }
    double s()    { return thermo.s;    }
    double e()    { return thermo.e;    }
    double rhoB() { return thermo.rhoB; }
    double rhoS() { return thermo.rhoS; }
    double rhoQ() { return thermo.rhoQ; }
    double w()    { return thermo.w;    }
    double A()    { return thermo.A;    }
    double cs2()  { return thermo.cs2;  }
    double dwds() { return thermo.dwds; }
    double dwdB() { return thermo.dwdB; }
    double dwdS() { return thermo.dwdS; }
    double dwdQ() { return thermo.dwdQ; }

    string get_current_eos_name() { return thermo.eos_name; }

};
}

/// @brief Overload the << operator for the Particle class
/// @details This allows for printing the contents of a Particle object to an
/// ostream object (typically for streaming to standard output or a file).
/// @tparam D The dimensionality of the simulation
/// @param os ostream object
/// @param p particle object
/// @return ostream object
template<unsigned int D>
inline ostream& operator<<( ostream& os, const ccake::Particle<D>& p ){
  os << "ID.............: " << p.ID << endl;
  os << "r..............: " << p.r << endl;
  os << "input.s........: " << p.input.s << endl;
  os << "input.e........: " << p.input.e << endl;
  os << "input.rhoB.....: " << p.input.rhoB << endl;
  os << "input.rhoQ.....: " << p.input.rhoQ << endl;
  os << "input.rhoS.....: " << p.input.rhoS << endl;
  os << "smoothed.s.....: " << p.smoothed.s << endl;
  os << "smoothed.e.....: " << p.smoothed.e << endl;
  os << "smoothed.rhoB..: " << p.smoothed.rhoB << endl;
  os << "smoothed.rhoQ..: " << p.smoothed.rhoQ << endl;
  os << "smoothed.rhoS..: " << p.smoothed.rhoS << endl;
  os << "specific.s.....: " << p.specific.s << endl;
  os << "specific.e.....: " << p.specific.e << endl;
  os << "specific.rhoB..: " << p.specific.rhoB << endl;
  os << "specific.rhoQ..: " << p.specific.rhoQ << endl;
  os << "specific.rhoS..: " << p.specific.rhoS << endl;
  os << "d_dt_spec.s....: " << p.d_dt_spec.s << endl;
  os << "d_dt_spec.e....: " << p.d_dt_spec.e << endl;
  os << "d_dt_spec.rhoB.: " << p.d_dt_spec.rhoB << endl;
  os << "d_dt_spec.rhoQ.: " << p.d_dt_spec.rhoQ << endl;
  os << "d_dt_spec.rhoS.: " << p.d_dt_spec.rhoS << endl;
  os << "norm_spec.s....: " << p.norm_spec.s << endl;
  os << "norm_spec.e....: " << p.norm_spec.e << endl;
  os << "norm_spec.rhoB.: " << p.norm_spec.rhoB << endl;
  os << "norm_spec.rhoQ.: " << p.norm_spec.rhoQ << endl;
  os << "norm_spec.rhoS.: " << p.norm_spec.rhoS << endl;
  os << "thermo.T.......: " << p.thermo.T << endl;
  os << "thermo.muB.....: " << p.thermo.muB << endl;
  os << "thermo.muQ.....: " << p.thermo.muQ << endl;
  os << "thermo.muS.....: " << p.thermo.muS << endl;
  os << "thermo.e.......: " << p.thermo.e << endl;
  os << "thermo.s.......: " << p.thermo.s << endl;
  os << "thermo.rhoB....: " << p.thermo.rhoB << endl;
  os << "thermo.rhoQ....: " << p.thermo.rhoQ << endl;
  os << "thermo.rhoS....: " << p.thermo.rhoS << endl;
  os << "thermo.p.......: " << p.thermo.p << endl;
  os << "thermo.cs2.....: " << p.thermo.cs2 << endl;
  os << "thermo.w.......: " << p.thermo.w << endl;
  os << "thermo.A.......: " << p.thermo.A << endl;
  os << "thermo.dwds....: " << p.thermo.dwds << endl;
  os << "thermo.dwdB....: " << p.thermo.dwdB << endl;
  os << "thermo.dwdQ....: " << p.thermo.dwdQ << endl;
  os << "thermo.dwdS....: " << p.thermo.dwdS << endl;
  os << "hydro.t........: " << p.hydro.t << endl;
  os << "hydro.gamma....: " << p.hydro.gamma << endl;
  os << "hydro.theta....: " << p.hydro.theta << endl;
  os << "hydro.v........: " << p.hydro.v << endl;
  os << "hydro.u........: " << p.hydro.u << endl;
  os << "hydro.rho_Q_ext: " << p.hydro.rho_Q_ext << endl;
  os << "hydro.rho_S_ext: " << p.hydro.rho_S_ext << endl;
  os << "hydro.rho_B_ext: " << p.hydro.rho_B_ext << endl;
  os << "hydro.gradP....: " << p.hydro.gradP << endl;
  os << "hydro.gradE....: " << p.hydro.gradE << endl;
  os << "hydro.gradBulk.: " << p.hydro.gradBulk << endl;
  os << "hydro.divshear.: " << p.hydro.divshear << endl;
  os << "hydro.gradshear: " << p.hydro.gradshear << endl;
  os << "hydro.M_big_bulk: " << p.hydro.M_big_bulk << endl;
  os << "hydro.M_shv_nabla_u: " << p.hydro.M_shv_nabla_u << endl;
  os << "hydro.M_big_entropy: " << p.hydro.M_big_entropy << endl;
  os << "hydro.R_big_entropy.: " << p.hydro.R_big_entropy << endl;
  os << "hydro.R_big_bulk...: " << p.hydro.R_big_bulk << endl;
  os << "hydro.F_u.....: " << p.hydro.F_u << endl;
  os << "hydro.M_u.....: " << p.hydro.M_u << endl;
  os << "hydro.F_big_N: " << p.hydro.F_big_N << endl;
  os << "hydro.F_0i_shear: " << p.hydro.F_0i_shear << endl;
  os << "hydro.M_0i_shear: " << p.hydro.M_0i_shear << endl;
  os << "hydro.F_big_entropy.....: " << p.hydro.F_big_entropy << endl;
  os << "hydro.shv_nabla_u: " << p.hydro.shv_nabla_u << endl;
  os << "hydro.bulk....: " << p.hydro.bulk << endl;
  os << "hydro.bigBulk.: " << p.hydro.bigBulk << endl;
  os << "hydro.tau_Pi..: " << p.hydro.tau_Pi << endl;
  os << "hydro.tau_pi..: " << p.hydro.tau_pi << endl;
  os << "hydro.delta_PiPi: " << p.hydro.delta_PiPi << endl;
  os << "hydro.delta_pipi: " << p.hydro.delta_pipi << endl;
  os << "hydro.lambda_Pipi: " << p.hydro.lambda_Pipi << endl;
  os << "hydro.varsigma.: " << p.hydro.varsigma << endl;
  os << "hydro.phi1....: " << p.hydro.phi1 << endl;
  os << "hydro.phi3....: " << p.hydro.phi3 << endl;
  os << "hydro.phi6....: " << p.hydro.phi6 << endl;
  os << "hydro.phi7....: " << p.hydro.phi7 << endl;
  os << "hydro.eta_pi..: " << p.hydro.eta_pi << endl;
  os << "hydro.tau_pipi: " << p.hydro.tau_pipi << endl;
  os << "hydro.zeta_Pi....: " << p.hydro.zeta_Pi << endl;
  os << "hydro.sigma_star: " << p.hydro.sigma_star << endl;
  os << "hydro.sigma...: " << p.hydro.sigma << endl;
  os << "hydro.dbigBulk_dt: " << p.hydro.dbigBulk_dt << endl;
  os << "hydro.F_big_bulk: " << p.hydro.F_big_bulk << endl;
  os << "hydro.F_shv_nabla_u: " << p.hydro.F_shv_nabla_u << endl;
  os << "hydro.F_big_entropy.....: " << p.hydro.F_big_entropy << endl;
  os << "hydro.F_big_N.: " << p.hydro.F_big_N << endl;
  os << "hydro.F_u.....: " << p.hydro.F_u << endl;
  os << "hydro.F_big_N.: " << p.hydro.F_big_N << endl;
  os << "hydro.du_dt...: " << p.hydro.du_dt << endl;
  os << "hydro.d_bigshv_dt.: " << p.hydro.d_bigshv_dt << endl;
  os << "hydro.shv.....: " << p.hydro.shv << endl;
  os << "hydro.bigshv..: " << p.hydro.bigshv << endl;
  os << "hydro.F_big_shear: " << p.hydro.F_big_shear << endl;
  os << "hydro.M_big_shear: " << p.hydro.M_big_shear << endl;
  os << "hydro.R_big_shear: " << p.hydro.R_big_shear << endl;
  return os;
}

#endif
