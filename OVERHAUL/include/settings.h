#ifndef SETTINGS_H
#define SETTINGS_H

#include <algorithm>
#include <string>
#include <vector>

using std::string;

class Settings
{
  public:

    // default settings
    bool using_Gubser            = false;
    bool using_Gubser_with_shear = false;
    bool using_shear             = false;

    bool printing_to_txt         = true;
    bool printing_to_HDF         = true;

    static constexpr int VERBOSE = 5;

    // maximum upper limit for t
    static constexpr double tend = 50.02;
    
    // simulation settings
    double t0      = 0.0; // initial timestep
    double dt      = 0.0; // size of timestep
    double stepx   = 0.0; // dx [fm]
    double stepy   = 0.0; // dy [fm]
    double h       = 0.0; // SPH kernel scale [fm]

    double Freeze_Out_Temperature = 0.0;  // freeze-out temperature (at zero mu)
//    double efcheck                = 0.0;  // freeze-out energy density
//    double sfcheck                = 0.0;  // freeze-out entropy density


    // quantities read in from InputParameters.inp file
    vector<string> headers;

    string IC_type         = "";  // specify initial condition type
    string IC_option       = "";  // specify option for given initial condition type
    string IC_file         = "";  // specify option for given initial condition type
    string EoS_type        = "";  // specify equation of state type
    string EoS_option      = "";  // specify specific option for EOS
    string eta             = "";  // specify the shear viscosity type to use
    string etaOption       = "";  // specify necessary eta options
    string shearRelax      = "";  // specify which shear relaxation to use
    string zeta            = "";  // specify the bulk viscosity type to use
    string zetaOption      = "";  // specify necessary zeta options
    string bulkRelax       = "";  // specify which bulk relaxation to use
    string Freeze_Out_Type = "";  // which freeze-out criterion to use

    // allows for explicitly printing extra information about specific particles
    vector<int> particles_to_print;

    // make sure that all chosen settings make reasonable sense
    void check_consistency()
    {
      // put any necessary consistency checks here

      return;
    }

};

typedef std::shared_ptr<Settings> pSettings;  // smart pointer to settings object

#endif