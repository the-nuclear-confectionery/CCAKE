#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include <cmath>
#include "vector.h"
#include "matrix.h"
#include "LinkList.h"
#include "tables.h"
#include "particle.h"
#include "eos.h"

namespace rk
{
	void bsq_second_order( double dx, void (*derivatives)( LinkList<D> &linklist ),
                         LinkList<D> &linklist );
}

#endif
