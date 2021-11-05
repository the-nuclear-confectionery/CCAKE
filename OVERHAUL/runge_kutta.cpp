#include <cmath>
#include "vector.h"
#include "matrix.h"
#include "linklist.h"
#include "particle.h"
//#include "eos.h"
#include "new_eos/eos.h"
#include "equations_of_motion.h"
#include "system_state.h"
#include "runge_kutta.h"

namespace RK
{
  void bsq_second_order( double dx, EquationsOfMotion & eom,
                         SystemState & system, SPHWorkstation & ws )
  {
    system.rk2 = 1;
    double t0  = system.t;
    double E0  = system.Ez;

    // initialize quantities at current time step
    system.set_current_timestep_quantities();
  
    ////////////////////////////////////////////
    //    first step
    ////////////////////////////////////////////

    // compute derivatives
    eom.BSQshear(system, ws);

    // update quantities
    system.get_derivative_halfstep(dx);

    system.Ez = E0 + 0.5*dx*system.dEz;
    system.t  = t0 + 0.5*dx;

    ////////////////////////////////////////////
    //    second step
    ////////////////////////////////////////////

    // compute derivatives
    eom.BSQshear(system, ws);

    // update quantities
    system.get_derivative_fullstep(dx);

    system.Ez = E0 + dx*system.dEz;
    system.t  = t0 + dx;

    return;
	}

}