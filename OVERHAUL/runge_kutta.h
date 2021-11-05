#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include <cmath>
#include "vector.h"
#include "matrix.h"
#include "linklist.h"
#include "particle.h"
//#include "eos.h"
#include "new_eos/eos.h"
#include "equations_of_motion.h"
#include "system_state.h"

namespace RK
{
	void bsq_second_order( double dx, EquationsOfMotion & eom, SystemState & system,
                          SPHWorkstation & ws );
}

#endif
