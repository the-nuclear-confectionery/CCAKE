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


    // quantities read in from InputParameters.inp file
    vector<string> headers;

    //------------------------------------
    // initial conditions
    string IC_type                   = "";
    string IC_option                 = "";
    string IC_file                   = "";

    //------------------------------------
    // equation of state
    string EoS_type                  = "";
    string EoS_option                = "";

    //------------------------------------
    // transport coefficients
    //  - shear quantities
    string etaMode                   = "";
    double constant_eta_over_s       = 0.0;
    string shearRelaxMode            = "";

    //  - bulk quantities
    string zetaMode                  = "";
    double constant_zeta_over_s      = 0.0;
    double cs2_dependent_zeta_A      = 0.0;
    double cs2_dependent_zeta_p      = 0.0;
    string bulkRelaxMode             = "";

    //------------------------------------
    // freeze out
    string Freeze_Out_Type           = "";



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