#ifndef SYSTEM_STATE_H
#define SYSTEM_STATE_H

#include <string>
#include <vector>

// functions calls to static EoS C library
#include <lib.h>
#include "eos_delaunay/eos_delaunay.h"

#include "vector.h"
#include "eos.h"
#include "matrix.h"
#include "particle.h"
#include "linklist.h"
#include "kernel.h"
#include "settings.h"

using std::string;
using std::vector;

class SystemState
{
friend class InputOutput;

public:

  void set_equation_of_state( EquationOfState & eos );
  void set_settings( Settings & settings );

  void compute_Knudsen_numbers();
  void compute_inverse_Reynolds_numbers();
  void compute_spatial_eccentricities( const vector<int> & orders );
  void compute_momentum_eccentricities( const vector<int> & orders );
  //double entropytotal;
  double S, S0;
  double t, dt;
  double Btotal, Btotal0;
  double Stotal, Stotal0;
  double Qtotal, Qtotal0;
  double E, Ez, E0, Eloss, dEz, Etot;
  int number_part, _n, N, rk2;
  double efcheck, sfcheck, freezeoutT;

  int etaconst, visc, cfon;
  double bvf, svf, zTc, sTc, zwidth;

  //void check_BSQ_E_conservation();
  //void check_BSQ_charge_conservation();

private:

  LinkList linklist;  // can we please name it something else?
  // vector<vector<int> > neighbors;  //?

  //EquationOfState eos;
  EquationOfState * eosPtr = nullptr;

  Settings * settingsPtr = nullptr;

  vector< Particle > particles;

  // creating vectors of vectors of the derivatives at each step
  vector<double> etasigma0;
  vector<double> Bulk0;

  vector< Vector<double,2> > u0;
  vector< Vector<double,2> > r0;

  vector< Matrix <double,2,2> > shv0;

  vector<int> list;


public:

  void initialize();
  void initialize_linklist(vector<Particle> & particles) { linklist.initialize(particles); }
  void BSQSimulation( double dt, LinkList & linklist );
  void BSQshear( LinkList & linklist );
  void check_BSQ_energy_conservation();
  void check_BSQ_entropy_conservation();
  void check_BSQ_charge_conservation();
  void smooth_fields(int a, bool init_mode = false);
  void smooth_gradients( int a, double tin, int & count );


  void bsqsvconservation();
  void conservation_entropy();
  void conservation_BSQ();
  void bsqsvconservation_E();
  void bsqsv_set();
  void bsqsvconservation_Ez();
  void setshear();

  void initialize_entropy_and_charge_densities();
  void initial_smoothing();

  int n(){ return _n; }


  // these routines are called in runge kutta
  void set_current_timestep_quantities();
  void get_derivative_halfstep(double dx);
  void get_derivative_fullstep(double dx);


};

#endif