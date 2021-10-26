#ifndef SYSTEM_STATE_H
#define SYSTEM_STATE_H

#include <string>
#include <vector>

#include "particle.h"
#include "runge_kutta.h"

using std::string;
using std::vector;

class SystemState
{
friend class InputOutput;

public:

  void compute_Knudsen_numbers();
  void compute_inverse_Reynolds_numbers();
  void compute_spatial_eccentricities( const vector<int> & orders );
  void compute_momentum_eccentricities( const vector<int> & orders );

  //void check_BSQ_E_conservation();
  //void check_BSQ_charge_conservation();

private:

  LinkList linklist;  // can we please name it something else?
  // vector<vector<int> > neighbors;  //?

  //EquationOfState eos;
  EquationOfState * eosPtr = nullptr;

  vector< Particle > particles;

  //RK runge_kutta_solver;

  //Evolve evolver;


public:

  void initialize();
  void BSQSimulation( double dt, LinkList & linklist );
  void BSQshear( LinkList & linklist );
  void check_BSQ_energy_conservation();
  void check_BSQ_entropy_conservation();
  void check_BSQ_charge_conservation();
  void smooth_fields(int a, bool init_mode = false);
  void smooth_gradients( int a, double tin, int & count );


  //void bsqsvconservation();
  //void conservation_entropy();
  //void conservation_BSQ();
  //void bsqsvconservation_E();
  void bsqsv_set();
  //void bsqsvconservation_Ez();
  void setshear();




};

#endif