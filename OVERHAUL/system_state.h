#include <string>
#include <vector>

#include "particle.h"

using std::string;
using std::vector;

class SystemState
{

public:

  void compute_Knudsen_numbers();
  void compute_inverse_Reynolds_numbers();
  void compute_spatial_eccentricities( const vector<int> & orders );
  void compute_momentum_eccentricities( const vector<int> & orders );

  void check_BSQ_E_conservation();
  void check_BSQ_charge_conservation();

private:

  LinkList linklist;  // can we please name it something else?
  // vector<vector<int> > neighbors;  //?

  vector< Particle > particles;

  rk::RK runge_kutta_solver;

  Evolve evolver;


}