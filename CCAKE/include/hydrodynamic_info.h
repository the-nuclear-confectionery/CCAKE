#ifndef HYDRODYNAMIC_INFO_H
#define HYDRODYNAMIC_INFO_H

#include <memory>

#include "matrix.h"
#include "vector.h"

/// \brief Strut to store hydrodynamic information
/// \tparam D Dimensionality of the problem
/// \details This struct is used to store hydrodynamic information for each
/// particle.
///
/// NOTE: Beware these information needs to be passed to the cabana AoSoA. To this
/// end, we define helper macros at the bottom of this file. If one wants to
/// add a new member to this struct, one needs to add it to the macros as well.
/// Also, remmember to pass new datamembers to cabana AoSoA in
/// `SystemState::allocate_cababa_particles()`
namespace ccake
{
template<unsigned int D>
struct hydrodynamic_info
{
  bool print_particle    = false;

  int ID                 = -1;  ///< for debugging purposes only
  double t               = 0.0; ///< current time in hydro simulation

  double Agam            = 0.0; ///< Agam = A ///TODO: check this
  double Agam2           = 0.0; ///< Agam2 = -Atot ///TODO: check this
  double shv33           = 0.0; ///< shear viscosity pi^33 (or pi_33) ///TODO: check this

  double gamma           = 0.0; ///< Lorentz factor
  double Bulk            = 0.0; ///< Bulk Viscosity weight ///TODO: what is this?
  double bigPI           = 0.0; ///< total bulk viscosity
  double C               = 0.0;
  double tauRelax        = 0.0; ///< Bulk Relaxation time
  double stauRelax       = 0.0; ///< Shear Relxation time
  double zeta            = 0.0; ///< bulk coefficient
  double setas           = 0.0; ///< half Shear coefficient eta ///TODO: check this
  double Ctot            = 0.0; ///< See Jaki's notes, Eq. (270)
  double Btot            = 0.0; ///< See Jaki's notes, Eq. (274)

  double sigma           = 0.0; ///< especific volume
  double dsigma_dt       = 0.0; ///< derivative of especific volume

  double g2              = 0.0; ///< gamma^2
  double g3              = 0.0; ///< gamma^3
  double gt              = 0.0; ///< gamma*tau
  double eta_o_tau       = 0.0; ///< shear visc coeff eta/tau/2 \\\TODO: check this
  double dwdsT1          = 0.0; ///< 1 -  (1/T) dw/ds
  double sigl            = 0.0; ///< (1/sigma^*) d sigma^*/dt - 1/tau

  double varsigma        = 0.0; ///< defined to be p + Pi (pressure + bulk)

  // vector members
  Vector<double,D> v;                     ///< velocity
  Vector<double,D> u;                     ///< relativistic velocity

  Vector<double,D> gradP;                 ///< Gradient of Pressure
  Vector<double,D> gradBulk;              ///< Gradient of Bulk Viscosity
  Vector<double,D> divshear, gradshear;

  Matrix<double,D,D> Imat;
  Matrix<double,D,D> gradV, gradU;        // Gradient of velocity needed for shear
  Matrix<double,D,D> uu, pimin, piu, piutot;
  Matrix<double,D+1,D+1> shv;
  Matrix<double,D,D> shv1, shv2, shv3, shv4;


  // quantities (possibly?) set in EoM classes
  double bigtheta        = 0.0;
  double inside          = 0.0;
  double div_u           = 0.0;


  // derivatives
  double dBulk_dt        = 0.0;

  Vector<double, D> du_dt;
  Matrix<double, D, D> dshv_dt;

};

namespace hydro_info
{

enum hydro_scalar_info
{
  t,
  Agam,
  Agam2,
  shv33,
  gamma,
  Bulk,
  bigPI,
  C,
  tauRelax,
  stauRelax,
  zeta,
  setas,
  Ctot,
  Btot,
  sigma,
  dsigma_dt,
  g2,
  g3,
  gt,
  eta_o_tau,
  dwdsT1,
  sigl,
  varsigma,
  bigtheta,
  inside,
  div_u,
  dBulk_dt,
  NUM_HYDRO_SCALAR_INFO
};
#define HYDRO_SCALAR_INFO double[ccake::hydro_info::NUM_HYDRO_SCALAR_INFO]
enum hydro_vector_info
{
  v,
  u,
  gradP,
  gradBulk,
  divshear,
  gradshear,
  du_dt,
  NUM_HYDRO_VECTOR_INFO
};
#define HYDRO_VECTOR_INFO double[ccake::hydro_info::NUM_HYDRO_VECTOR_INFO][3]
enum hydro_space_matrix_info
{
  Imat,
  gradV,
  gradU,
  uu,
  pimin,
  piu,
  piutot,
  shv1,
  shv2,
  shv3,
  shv4,
  dshv_dt,
  NUM_HYDRO_SPACE_MATRIX_INFO
};
#define HYDRO_SPACE_MATRIX_INFO double[ccake::hydro_info::NUM_HYDRO_SPACE_MATRIX_INFO][3][3]
enum hydro_spacetime_matrix_info
{
  shv,
  NUM_HYDRO_SPACETIME_MATRIX_INFO
};
#define HYDRO_SPACETIME_MATRIX_INFO double[ccake::hydro_info::NUM_HYDRO_SPACETIME_MATRIX_INFO][4][4]

}}

#endif