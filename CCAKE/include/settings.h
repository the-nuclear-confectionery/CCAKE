#ifndef SETTINGS_H
#define SETTINGS_H

#include <algorithm>
#include <string>
#include <vector>
#include <filesystem>

#include "constants.h"
#include "formatted_output.h"

using std::string;
namespace fs = std::filesystem;

class Settings
{
  public:

    //==========================================================================
    // default/global settings
    bool using_Gubser                 = false;  //TODO: Remove this
    bool using_Gubser_with_shear      = false;  //TODO: Remove this
    bool using_shear                  = false;  //TODO: Whatever to use shear or
                                                // not should be decided by the
                                                // eta mode, not here
    bool initializing_with_full_Tmunu = false;  // whether to initialize Pi from
                                                // varsigma - p or not
    //==========================================================================
    // I/O settings
    bool printing_to_txt              = true;
    bool printing_to_HDF              = true;

    static constexpr int VERBOSE      = 5;

    //==========================================================================
    // hydrodynamics settings

    // restrict number of timesteps, if desired (negative means no restriction)
    static constexpr int max_number_of_timesteps = -1;

    // maximum upper limit for t
    static constexpr double tend      = 50.02;

    // simulation settings
    double stepx                      = 0.0; // dx [fm]
    double stepy                      = 0.0; // dy [fm]
    double stepEta                    = 0.0; // d\eta
    double xmin                       = 0.0; // minimum x [fm]
    double ymin                       = 0.0; // minimum y [fm]
    double etamin                     = 0.0; // minimum eta

    vector<string> headers;
    // allows for explicitly printing extra information about specific particles
    vector<int> particles_to_print;

    //==========================================================================
    // input parameter settings
    // quantities read in from InputParameters.inp file


    //------------------------------------
    // initial conditions
    string IC_type              = "";    ///< Type of initial conditions
                                         ///  ("ccake", "iccing", "hdf5 )
                                         //TODO: Implement hdf5 reader
    string IC_option            = "";    ///< Deprecated. Replaced by the enabling
                                         ///  disabling B, S or Q charges in hydro conf
                                         ///  //TODO: Remove this
    fs::path IC_file            = "";    ///< Path to the initial conditions file
    double t0                   = 0.0;   ///< Proper time of the initial conditions
    unsigned int dim            = 2;     ///< Dimension of the initial conditions (1, 2, or 3)
    bool input_as_entropy       = false; ///< Whether the initial conditions are
                                         /// given as entropy density (true) or energy density

    //------------------------------------
    // simulations parameters
    double dt                   = 0.0;   ///< size of timestep
    double hT                   = 0.0;   ///< SPH kernel scale in transverse direction [fm]
    double hEta                 = 0.0;   ///< SPH kernel scale in longitudinal direction
    std::string kernel_type     = "";    ///< Which SPH kernel to use (cubic, quartic, quintic).
                                         ///  Currently only cubic is supported
                                         //TODO: Implement other kernels
    double e_cutoff              = 0.0;  ///< energy density below which
                                         ///particles are removed [GeV/fm^3]

    //------------------------------------
    // buffer particles settings
    bool buffer_event           = true;   ///< add a buffer around event to
                                          ///  stabilize evolution
    bool circular_buffer        = true;   ///< whether to buffer with entire
                                          ///  grid or just circular padding
    double padding_thickness    = 0.1;    ///< if circular_buffer == true,
                                          ///  buffer_radius specifies the
                                          ///  fractional amount of padding
                                          ///  to add beyond the point with
                                          ///  maximum distance from origin

    //------------------------------------
    // Particlization settings
    bool particlization_enabled = false;  ///< whether to use the "old" or "new"
                                          ///  particlization scheme
    string Freeze_Out_Type      = "";     ///< Type of freeze-out ("fixed_T") //TODO: Implement other types (e.g., fixed energy and fixed tau)
    double Freeze_Out_Temperature = 0.0;  ///< freeze-out temp. (at zero mu)

    //------------------------------------
    // equation of state
    string eos_type             = "";     ///< Type of equation of state ("conformal" or "table")
    string eos_path             = "";     ///< If "table", path to the equation of state file

    //------------------------------------
    // hydrodynamics configuration
    bool baryon_charge_enabled  = true;   ///< whether to include baryon charge
                                          ///  in the hydrodynamics
    bool strange_charge_enabled = true;   ///< whether to include strangess charge
                                          ///  in the hydrodynamics
    bool electric_charge_enabled = true;  ///< whether to include electric charge
                                          ///  in the hydrodynamics
    //------------------------------------
    // transport coefficients
    //  - shear quantities
    string etaMode              = "";
    double constant_eta_over_s  = 0.0;
    string shearRelaxMode       = "";

    //  - bulk quantities
    string zetaMode             = "";
    double constant_zeta_over_s = 0.0;
    double cs2_dependent_zeta_A = 0.0;
    double cs2_dependent_zeta_p = 0.0;
    string bulkRelaxMode        = "";
    bool modulate_zeta_with_tanh = true;   // forces zeta/s to decrease
                                           // smoothly to zero below
                                           // transition temperature

    //------------------------------------
    // Gubser settings
    string Gubser_BSQmode            = ""; //TODO: Remove this


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
        if ( eos_type != "conformal" )
        {
          std::cerr << "WARNING: Gubser initial conditions require a conformal "
                       "equation of state!  Switching to gas of massless gluons"
                       " and 2.5 massless quarks" << std::endl;
          eos_type = "conformal";
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