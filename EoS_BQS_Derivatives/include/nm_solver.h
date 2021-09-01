#ifndef NM_SOLVER_H
#define NM_SOLVER_H

#include "Functions.h"

//void solve( double & e0, double & B0, double & S0, double & Q0,
//			double & Tout, double & muBout, double & muSout, double & muQout );
void solve( double densities[], double sols[] );
void solve2( double densities[], double sols[], double minima[], double maxima[] );

#endif
