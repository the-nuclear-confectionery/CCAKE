#ifndef SETTINGS_H
#define SETTINGS_H

#include <algorithm>
#include <string>

using std::string;

class Settings
{
public:

  bool using_Gubser            = false;
  bool using_Gubser_with_shear = false;
  bool using_shear             = false;

  /* DOUBLES */
  // maximum upper limit for t
  static constexpr double tend=50.02;
  
  double step = 0.0;
  double t0 = 0.0;
  double t = 0.0;
  double dt = 0.0;
  double stepx = 0.0;
  double stepy = 0.0; // initial grid coordinate spacings
  double stepr = 0.0;
  double stepphi = 0.0;
  double _h = 0.0;
  double Freeze_Out_Temperature = 0.0;

  double gd2 = 0.0;
  double factor = 0.0;

  double efcheck = 0.0;
  double sfcheck = 0.0;

  double E1 = 0.0;
  double E2 = 0.0;

/********************************************************************************/


/* INTS */
  static constexpr int VERBOSE = 5;
  //static constexpr double q=1.;

  //int _n  = 0;
  int start = 0;
  int end   = 0;
  int fnum  = 0;
  int qmf   = 3;  // if==1 quantum mechanicanical corrections to the flow or added,
                  // if==0 no corrections are included
  //int number_part = 0;
  //int rk2 = 0;
  int gtyp = 0;
  //int cfon = 0;
  int visc = 0; // visc=0 for ideal
            // visc=1 for bulk,
            // visc=2 for shear
            // visc=3 for bulk+shear
            // visc=4 for bulk+shear+BSQ

  int steps    = 0;
  //int first  = 0;
  int average  = 0;
  int lowT     = 0;
  int etaconst = 0;

  int fcount = 0;
  int cevent = 0;
  /********************************************************************************/


  /* STRINGS */

  string eos_s = "";
  string eos_p = "";
  string ebe_folder = "";
  string eost = "";
  string IC_type = ""; // specify initial condition type
  string IC_option = ""; // specify option for given initial condition type
  string IC_file = ""; // specify option for given initial condition type
  string EoS_type = ""; // specify equation of state type
  string EoS_option = ""; // specify specifc option for EOS
  // there should an associated EoS directory with tables
  string eta = ""; // specificy the shear viscosity type to use
  string etaOption = ""; // specify necessary eta options
  string shearRelax = ""; // specify which shear relaxation to use
  string zeta = ""; // specificy the bulk viscosity type to use
  string zetaOption = ""; // specify necessary zeta options
  string bulkRelax = ""; // specify which bulk relaxation to use
  string Freeze_Out_Type = "";
  string initial_coordinate_distribution = "";


/********************************************************************************/


/* VECTOR of headers */
  vector<string> headers;
/********************************************************************************/


  // allows for explicitly printing tons of information about specific particles
  vector<int> particles_to_print;


  // make sure that all chosen settings make reasonable sense
  void check_consistency()
  {
    // put any necessary consistency checks here

    return;
  }

  inline bool print_particle( int i )
  {
    return std::binary_search( particles_to_print.begin(),
                               particles_to_print.end(), i );
  }

};





#endif