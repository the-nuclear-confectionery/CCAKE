#ifndef HYDRODYNAMIC_INFO_H
#define HYDRODYNAMIC_INFO_H

#include <memory>

#include "matrix.h"
#include "vector.h"

// hydrodynamic quantities (with default initialization)
struct hydrodynamic_info
{
  bool print_particle    = false;

  int ID                 = -1;  // for debugging purposes only
  double t               = 0.0; // current time in hydro simulation

  double Agam            = 0.0;
  double Agam2           = 0.0;
  double shv33           = 0.0;

  double gamma           = 0.0; // Lorentz factor
  double Bulk            = 0.0; // Bulk Viscosity weight
  double bigPI           = 0.0; // total bulk viscosity
  double C               = 0.0;
  double tauRelax        = 0.0; // Bulk Relaxation time
  double stauRelax       = 0.0; // Shear Relxation time
  double zeta            = 0.0; // bulk coefficient
  double setas           = 0.0; 
  double Ctot            = 0.0;
  double Btot            = 0.0;

  double sigma           = 0.0; // especific volume
  double dsigma_dt       = 0.0; // derivative of especific volume

  double g2              = 0.0;
  double g3              = 0.0;
  double gt              = 0.0;
  double eta_o_tau       = 0.0;
  double dwdsT1          = 0.0;
  double sigl            = 0.0;

  double varsigma        = 0.0; // defined to be p + Pi (pressure + bulk)

  // vector members
  Vector<double,2> v;                     // velocity
  Vector<double,2> u;                     // relativistic velocity

  Vector<double,2> gradP;                 // Gradient of Pressure
  Vector<double,2> gradBulk;              // Gradient of Bulk Viscosity
  Vector<double,2> divshear, gradshear;


  // matrix members
  Matrix<double,2,2> Imat;
  Matrix<double,2,2> gradV, gradU;        // Gradient of velocity needed for shear
  Matrix<double,2,2> uu, pimin, piu, piutot;
  Matrix<double,3,3> shv;
  Matrix<double,2,2> shv1, shv2, shv3, shv4;


  // quantities (possibly?) set in EoM classes
  double bigtheta        = 0.0;
  double inside          = 0.0;
  double div_u           = 0.0;


  // derivatives
  double dBulk_dt        = 0.0;

  Vector<double, 2> du_dt;
  Matrix<double, 2, 2> dshv_dt;


};

typedef std::shared_ptr<hydrodynamic_info> pHydrodynamic_info;


#endif