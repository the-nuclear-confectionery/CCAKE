#ifndef EOS_TABLE_H
#define EOS_TABLE_H

#include <string>

// functions calls to static EoS C library
#include <lib.h>
#include "eos_delaunay/eos_delaunay.h"
#include "interpolatorND/interpolatorND.h"
#include "constants.h"
#include "eos_base.h"
#include "eos_header.h"

class EoS_table: public EoS_base
{
private:
//	const double hc = constants::hbarc_MeVfm;
//
//  double c, T0, muB0, muS0, muQ0;

  static InterpolatorND<4> equation_of_state_table;

  void init_grid_ranges_only(string quantityFile, string derivFile);


  void get_eBSQ_densities_from_interpolator( double point[], double densities[] );
  void get_sBSQ_densities_from_interpolator( double point[], double densities[] );

public:
  // default constructor/destructor
  EoS_table(){}
  ~EoS_table(){}

  // constructor from files containing tables
  EoS_table( string quantityFile, string derivFile );

  void get_eBSQ( double point[], double results[] );
  void get_sBSQ( double point[], double results[] );
  void get_full_thermo( double point[], double results[] );
};

#endif