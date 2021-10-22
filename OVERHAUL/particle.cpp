#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "particle.h"

// Default constructor
particle::particle()
{
  r.resize(2, 0.0);
  v.resize(2, 0.0);
  u.resize(2, 0.0);
}

void particle::set_equation_of_state( EquationOfState & eos_in )
{
  eos = eos_in;
}


// 
double particle::locate_phase_diagram_point_eBSQ(
                 double e_In, double rhoB_In, double rhoS_In, double rhoQ_In )
{
  double sVal = EOS.locate_phase_diagram_point_eBSQ(
                    e_In, rhoB_In, rhoS_In, rhoQ_In );
  thermo.set(EOS);
  return sVal;
}


double particle::locate_phase_diagram_point_eBSQ(double e_In)
{
	return EOSs_out( e_In, 0.0, 0.0, 0.0 );
}

void particle::locate_phase_diagram_point_sBSQ(
                 double s_In, double rhoB_In, double rhoS_In, double rhoQ_In )
{
  bool update_s_success
       = EOS.locate_phase_diagram_point_sBSQ( s_In, rhoB_In, rhoS_In, rhoQ_In );
  thermo.set(EOS);

  return;
}

void particle::locate_phase_diagram_point_sBSQ(double s_In)
{
  locate_phase_diagram_point_sBSQ(s_In, 0.0, 0.0, 0.0 );
  return;
}