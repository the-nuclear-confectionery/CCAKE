#include <cmath>
#include "vector.h"
#include "matrix.h"
#include "linklist.h"
#include "particle.h"
#include "eos.h"
#include "equations_of_motion.h"
#include "system_state.h"
#include "runge_kutta.h"

namespace RK
{
  void bsq_second_order( double dx, EquationsOfMotion & eom, SystemState & system,
                          SPHWorkstation & ws )
  {
//    int N = system.n();
    N = system.n();


    double E0, t0;

    system.rk2 = 1;
    t0           = system.t;

    // initialize quantities at current time step
    system.set_current_timestep_quantities();
  

    E0 = system.Ez;

    ////////////////////////////////////////////
    //    first step
    ////////////////////////////////////////////

    // compute derivatives
    //(*derivatives)(linklist);
    eom.BSQshear(system, ws);

    // update quantities
    system.get_derivative_halfstep(dx);

    system.Ez = E0 + 0.5*dx*system.dEz;
    system.t  = t0 + 0.5*dx;

    ////////////////////////////////////////////
    //    second step
    ////////////////////////////////////////////

    // compute derivatives
    //(*derivatives)(linklist);
    eom.BSQshear(system, ws);

    // update quantities
    system.get_derivative_fullstep(dx);

    system.Ez = E0 + dx*system.dEz;
    system.t  = t0 + dx;

	}

}