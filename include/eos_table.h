#ifndef EOS_TABLE_H
#define EOS_TABLE_H

#include <algorithm>
#include <iterator>
#include <string>

#include "constants.h"
#include "eos_base.h"
#include "eos_gap.h"
#include "eos_header.h"
#include "interpolatorND.h"

using namespace constants;

class EoS_table: public EoS_base
{
private:
  bool debug_mode = false;
  bool using_HDF  = false;

  string equation_of_state_table_filename = "";

  std::shared_ptr<InterpolatorND<4>> equation_of_state_table;

  void get_eBSQ_densities_from_interpolator( double point[], double densities[] );
  void get_sBSQ_densities_from_interpolator( double point[], double densities[] );

  void get_eBSQ_safe( const double point[], double results[] );
  void get_sBSQ_safe( const double point[], double results[] );
  void get_full_thermo_safe( const double point[], double results[] );

  bool point_is_in_range_no_ext( double setT, double setmuB, double setmuQ, double setmuS );


public:
  GapTable gap_table; ///< auxiliary scattered-points table for the spinodal gap region

  // default constructor/destructor
  EoS_table(){}
  virtual ~EoS_table(){}

  // constructor from files containing tables
  EoS_table( string eos_path, bool normalize_by_T = true );

  /// true when the table carries the optional 18th field chiBB_nc (the non-critical
  /// chi_BB baseline used as a LOCAL xi_0 by the critical bulk scaling).
  bool has_chiBB_nc = false;

  void get_eBSQ( double point[], double results[] );
  void get_sBSQ( double point[], double results[] );
  void get_full_thermo( double point[], double results[] );

  void set_debug_mode( bool debug_mode_in ) { debug_mode = debug_mode_in; }
};

#endif