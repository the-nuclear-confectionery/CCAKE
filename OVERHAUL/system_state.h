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
friend class EquationsOfMotion;
friend class SPHWorkstation;

public:

  SystemState(){}
  ~SystemState(){}

  void set_EquationOfStatePtr( EquationOfState * eosPtr_in );
  void set_SettingsPtr(Settings * settingsPtr_in);

  //void compute_Knudsen_numbers();
  //void compute_inverse_Reynolds_numbers();
  //void compute_spatial_eccentricities( const vector<int> & orders );
  //void compute_momentum_eccentricities( const vector<int> & orders );
  //double entropytotal;


  /* DOUBLES */
  double _h      = 0.0;
  double S       = 0.0;
  double S0      = 0.0;
  double t       = 0.0;
  double dt      = 0.0;
  double Btotal  = 0.0;
  double Btotal0 = 0.0;
  double Stotal  = 0.0;
  double Stotal0 = 0.0;
  double Qtotal  = 0.0;
  double Qtotal0 = 0.0;

  double E          = 0.0;
  double Ez         = 0.0;
  double E0         = 0.0;
  double Eloss      = 0.0;
  double dEz        = 0.0;
  double Etot       = 0.0;
  double efcheck    = 0.0;
  double sfcheck    = 0.0;
  double freezeoutT = 0.0;
  double bvf        = 0.0;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  double svf        = 0.04; // HARDCODED FOR NOW
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  double zTc        = 0.0;
  double sTc        = 0.0;
  double zwidth     = 0.0;

 /* INTS */
  int number_part = 0;
  int _n          = 0;
  int N           = 0;
  int rk2         = 0;

  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  int etaconst = 1; // HARDCODED FOR NOW
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  int visc     = 4; // HARDCODED FOR NOW
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  int cfon     = 1; // HARDCODED FOR NOW
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


private:
/* POINTERS */
  EquationOfState * eosPtr = nullptr;
  Settings * settingsPtr   = nullptr;
/**********************************************************/



  LinkList linklist;  // can we please name it something else?
  // vector<vector<int> > neighbors;  //?

  //EquationOfState eos;
  vector<int> particles_out_of_grid;
  vector< Particle > particles;

  // creating vectors of vectors of the derivatives at each step
  vector<double> etasigma0;
  vector<double> Bulk0;

  vector< Vector<double,2> > u0;
  vector< Vector<double,2> > r0;

  vector< Matrix <double,2,2> > shv0;

  vector<int> list;

  // freeze-out related quantities
  int cf            = 0;
  int frzc          = 0;
  double avgetasig  = 0;
  double cs2        = 0;
  double tau        = 0;
  double taup       = 0;
  double taupp      = 0;
  double wfz        = 0;


  vector<double> divTtemp;
  vector<double> gsub;
  vector<double> bulksub;
  vector<double> swsub;
  vector<double> shear33sub;
  vector<double> tlist;
  vector<double> sFO;         //entropy at freezeout
  vector<double> Tfluc;
  vector< Vector<double,2> > divT;
  vector< Vector<double,2> > rsub;
  vector< Vector<double,2> > uout;

  //matrix pointer
  vector< Matrix<double,3,3> > shearsub;
  



public:

  void initialize();
  void reset_linklist() { linklist.reset(); }
  void initialize_linklist();
  void check_BSQ_energy_conservation();
  void check_BSQ_entropy_conservation();
  void check_BSQ_charge_conservation();

  void bsqsvconservation();
  void conservation_entropy();
  void conservation_BSQ();
  void bsqsvconservation_E();
  void bsqsv_set();
  void bsqsvconservation_Ez();

  int n(){ return _n; }


  // these routines are called in runge kutta
  void set_current_timestep_quantities();
  void get_derivative_step( double dx );

  // put FO routines in SystemState for time being
  void bsqsvfreezeout( int curfrz );
  void bsqsvinterpolate( int curfrz );

};

#endif