#ifndef SETTINGS_H
#define SETTINGS_H

#include <algorithm>
#include <string>
#include <vector>

#include "constants.h"
#include "formatted_output.h"

using std::string;

class Settings
{
  public:

    //==========================================================================
    // default/global settings
    bool using_Gubser                 = false;
    bool using_Gubser_with_shear      = false;
    bool using_shear                  = false;
    bool initializing_with_full_Tmunu = false;  // whether to initialize Pi from
                                                // varsigma - p or not

    //==========================================================================
    // event buffer settings
    bool buffer_event                 = false;   // add a buffer around event to
                                                // stabilize evolution
    bool circular_buffer              = true;   // whether to buffer with entire
                                                // grid or just circular padding
    double padding_thickness          = 0.1;    // if circular_buffer == true,
                                                // buffer_radius specifies the
                                                // fractional amount of padding
                                                // to add beyond the point with
                                                // maximum distance from origin
    

    //==========================================================================
    // transport coefficients settings
    bool modulate_zeta_with_tanh      = true;   // forces zeta/s to decrease
                                                // smoothly to zero below
                                                // transition temperature

    //==========================================================================
    // I/O settings
    bool printing_to_txt              = true;
    bool printing_to_HDF              = true;

    static constexpr int VERBOSE      = 0;

    //==========================================================================
    // hydrodynamics settings

    // restrict number of timesteps, if desired (negative means no restriction)
    static constexpr int max_number_of_timesteps = -1;

    // maximum upper limit for t
    static constexpr double tend      = 50.02;
    
    // simulation settings
    double t0                         = 0.0; // initial timestep
    double dt                         = 0.0; // size of timestep
    double stepx                      = 0.0; // dx [fm]
    double stepy                      = 0.0; // dy [fm]
    double stepEta                    = 0.0; // d\eta
    double xmin                       = 0.0; // minimum x [fm]
    double ymin                       = 0.0; // minimum y [fm]
    double etamin                     = 0.0; // minimum eta
    double h                          = 0.0; // SPH kernel scale [fm]

    double e_cutoff                   = 0.0; // energy density below which
                                             // particles are removed [GeV/fm^3]
    double Freeze_Out_Temperature     = 0.0; // freeze-out temp. (at zero mu)


    //==========================================================================
    // input parameter settings

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
    string EoS_path                  = "";

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

    //------------------------------------
    // Gubser settings
    string Gubser_BSQmode            = "";


    // allows for explicitly printing extra information about specific particles
    vector<int> particles_to_print;


    // make sure that all chosen settings make reasonable sense
    void check_consistency()
    {
      formatted_output::update("Impose consistency checks");

      //------------------------------------
      // enforce appropriate settings for Gubser
      if ( IC_type == "Gubser" || IC_type == "Gubser_with_shear" )
      {
        using_Gubser = true;
        if ( IC_type == "Gubser_with_shear" )
          using_Gubser_with_shear = true;

        //------------------------------------
        // put a warning check here
        if ( EoS_type != "conformal" )
        {
          std::cerr << "WARNING: Gubser initial conditions require a conformal "
                       "equation of state!  Switching to gas of massless gluons"
                       " and 2.5 massless quarks" << std::endl;
          EoS_type = "conformal";
        }

        //------------------------------------
        // run Gubser indefinitely
        Freeze_Out_Temperature = 1e-10/constants::hbarc_MeVfm;

        //------------------------------------
        // Gubser shear viscosity settings
        etaMode = "constant";
        if ( IC_type == "Gubser" )
          constant_eta_over_s = 0.0;
        else if ( IC_type == "Gubser_with_shear" )
          constant_eta_over_s = 0.20;

        //------------------------------------
        // Gubser bulk viscosity settings
        zetaMode = "constant";
        constant_zeta_over_s = 0.0;
      }
      else if ( IC_type == "TECHQM" )
      {
        t0 = 0.6;  //fm/c
      }

      // if eta/s == 0 identically, set using_shear to false
      if ( etaMode == "constant" && constant_eta_over_s < 1e-10 )
        using_shear  = false;
      else
        using_shear  = true;

      return;
    }

};

typedef std::shared_ptr<Settings> pSettings;  // smart pointer to settings object

#endif