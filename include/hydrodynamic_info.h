#ifndef HYDRODYNAMIC_INFO_H
#define HYDRODYNAMIC_INFO_H

#include <memory>

#include "matrix.h"
#include "vector.h"


namespace ccake
{
/// @brief Strut to store hydrodynamic information
/// @tparam D Dimensionality of the problem
/// @details This struct is used to store hydrodynamic information for each
/// particle. How these quantities are calculated is up to the user. For a
/// typical reference, see the documentation of the EoM_default class.
/// @see EoM_default
/// @note Beware these information needs to be passed to the cabana AoSoA. To this
/// end, we define helper macros at the bottom of this file. If one wants to
/// add a new member to this struct, one needs to add it to the macros as well.
/// Also, remmember to pass new datamembers to cabana AoSoA in
/// `SystemState::allocate_cababa_particles()`
template<unsigned int D>
struct hydrodynamic_info
{
  bool print_particle    = false;

  int ID                 = -1;  ///< for debugging purposes only
  double t               = 0.0; ///< current time in hydro simulation
  double a               = 0.0; ///< 0 for IR, 1 for DNMR
  double gamma           = 0.0; ///< Lorentz factor
  double theta           = 0.0; ///< expansion rate nabla_mu u^mu
  double bulk            = 0.0; ///< total bulk viscosity
  double extensive_bulk         = 0.0; ///< extensive bulk viscosity
  double tau_Pi          = 0.0; ///< Bulk Relaxation time
  double tau_pi          = 0.0; ///< Shear Relxation time
  double delta_PiPi      = 0.0; ///< bulk delta transport coefficient
  double delta_pipi      = 0.0; ///< shear delta transport coefficient
  double lambda_Pipi     = 0.0; ///< bulk-shear lambda transport coefficient
  double lambda_piPi     = 0.0; ///< shear-bulk lambda transport coefficient
  double tau_pipi        = 0.0; ///< shear relaxation time
  double phi1            = 0.0; ///< bulk phi1 transport coefficient
  double phi3            = 0.0; ///< bulk phi3 transport coefficient
  double phi6            = 0.0; ///< shear phi6 transport coefficient
  double phi7            = 0.0; ///< shear phi7 transport coefficient
  double zeta_Pi         = 0.0; ///< bulk coefficient
  double eta_pi          = 0.0; ///< shear coefficient
  double sigma_lab      = 0.0; ///< specific volume in computational frame
  double sigma           = 0.0; ///< specific volume
  double shv_nabla_u     = 0.0; ///< pi_mu_nu nabla^mu u^nu = pi^mu_nu sigma^mu_nu

  double varsigma        = 0.0; ///< defined to be p + Pi (pressure + bulk) 
  double rho_Q_ext       = 0.0; ///< external charge density
  double rho_S_ext       = 0.0; ///< external strangeness density
  double rho_B_ext       = 0.0; ///< external baryon density


  // vector members
  Vector<double,D> v;                     ///< velocity
  Vector<double,D> u;                     ///< relativistic velocity
  Vector<double,D> j_ext;                 ///< external current

  Vector<double,D> gradP;                 ///< Gradient of Pressure
  Vector<double,D> gradE;                 ///< Gradient of Energy
  Vector<double,D> gradBulk;              ///< Gradient of Bulk Viscosity
  Vector<double,D> divshear, gradshear;


  Matrix<double,D,D> gradV, gradU;        // Gradient of velocity needed for shear
  Matrix<double,4,4> shv;

  // Auxiliary M matrices
  Vector<double,D> M_extensive_bulk;
  Vector<double,D> M_shv_nabla_u;
  Vector<double,D> M_extensive_entropy;
  Matrix<double,D,D> M_u;
  //one for each conserved charge (baryon, strangeness, electric charge)
  Matrix<double,D,3> M_extensive_N;
  //one for each independent component of the shear tensor, linearized index
  Matrix<double,2,3*D> M_extensive_shear;
  Matrix<double,2,3*3> R_extensive_shear;

  // Auxiliary R matrices
  Vector<double,3> R_extensive_entropy;
  Vector<double,3> R_extensive_bulk;
  Matrix<double,3,3> R_extensive_N;
  Matrix<double,D,3> R_0i_shear;
  Matrix<double,D,D> M_0i_shear;
  Matrix<double,D,3> R_u;

  // auxiliary F vectors
  double F_extensive_bulk;
  double F_shv_nabla_u;
  double F_extensive_entropy;
  Vector<double,D> F_u;
  Vector<double,D> F_0i_shear;
  //one for each conserved charge (baryon, strangeness, electric charge)
  Vector<double,3> F_extensive_N;
  //one for each independent component of the shear tensor
  Matrix<double,2,3> F_extensive_shear;
  Matrix<double,2,3> extensive_shv;


  // derivatives
  double d_extensive_bulk_dt        = 0.0;
  Vector<double, D> du_dt;
  Matrix<double, 2, 3> d_extensive_shv_dt;

};

namespace hydro_info
{

enum hydro_scalar_info
{
  t,
  bulk,
  extensive_bulk,
  a,
  rho_Q_ext,
  rho_S_ext,
  rho_B_ext,
  tau_Pi,
  tau_pi,
  zeta_Pi,
  sigma_lab,
  sigma,
  gamma,
  theta,
  delta_PiPi,
  delta_pipi,
  lambda_Pipi,
  phi1,
  phi3,
  phi6,
  phi7,
  eta_pi,
  tau_pipi,
  lambda_piPi,
  varsigma,
  div_u,
  d_extensive_bulk_dt,
  F_extensive_bulk,
  F_shv_nabla_u,
  F_extensive_entropy,
  shv_nabla_u,
  NUM_HYDRO_SCALAR_INFO
};
#define HYDRO_SCALAR_INFO double[ccake::hydro_info::NUM_HYDRO_SCALAR_INFO]
enum hydro_vector_info
{
  v,
  u,
  gradP,
  gradE,
  gradBulk,
  divshear,
  gradshear,
  M_extensive_bulk,
  M_shv_nabla_u,
  M_extensive_entropy,
  R_extensive_entropy,
  R_extensive_bulk,
  F_u,
  F_extensive_N,
  F_0i_shear,
  j_ext,
  du_dt,
  NUM_HYDRO_VECTOR_INFO
};
#define HYDRO_VECTOR_INFO double[ccake::hydro_info::NUM_HYDRO_VECTOR_INFO][3]
enum hydro_space_matrix_info
{
  gradV,
  M_u,
  M_extensive_N,
  R_extensive_N,
  R_u,
  R_0i_shear,
  M_0i_shear,
  NUM_HYDRO_SPACE_MATRIX_INFO
};
#define HYDRO_SPACE_MATRIX_INFO double[ccake::hydro_info::NUM_HYDRO_SPACE_MATRIX_INFO][3][3]
enum hydro_spacetime_matrix_info
{
  shv,
  NUM_HYDRO_SPACETIME_MATRIX_INFO
};
#define HYDRO_SPACETIME_MATRIX_INFO double[ccake::hydro_info::NUM_HYDRO_SPACETIME_MATRIX_INFO][4][4]
enum hydro_shear_aux_vector_info
{
  F_extensive_shear,
  extensive_shv,
  d_extensive_shv_dt,
  NUM_HYDRO_SHEAR_AUX_VECTOR_INFO
};
#define HYDRO_SHEAR_AUX_VECTOR_INFO double[ccake::hydro_info::NUM_HYDRO_SHEAR_AUX_VECTOR_INFO][2][3]
enum hydro_shear_aux_matrix_info
{
  M_extensive_shear,
  R_extensive_shear,
  NUM_HYDRO_SHEAR_AUX_MATRIX_INFO
};
#define HYDRO_SHEAR_AUX_MATRIX_INFO double[ccake::hydro_info::NUM_HYDRO_SHEAR_AUX_MATRIX_INFO][2][9]


}}

#endif