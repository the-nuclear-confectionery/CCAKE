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

/// @class Settings
/// @brief A class to hold the settings of the simulation.
/// @details This class holds the settings of the simulation, such as the
/// initial conditions, the equation of state, the transport coefficients etc
/// The settings are read from a YAML file and stored in this class.
/// @todo Many options are hardcoded and should be read from the YAML file.
/// Others are not used anymore. Others yet are only set, but functionality
/// never implemented. A major cleanup is needed.
class Settings
{
  public:

    //==========================================================================
    // default/global settings
    bool using_shear                  = false;  //TODO: Whatever to use shear or
                                                // not should be decided by the
                                                // eta mode, not here
    bool using_diffusion              = false;  //TODO: Whatever to use diffusion or
                                                // not should be decided by the
                                                // eta mode, not here
    bool initializing_with_full_Tmunu = false;  // whether to initialize Pi from
                                                // tmunu_trace - p or not
    //==========================================================================
    // I/O settings
    bool txt_evolution              = false;
    bool hdf_evolution              = false;
    bool print_conservation_status  = true;
    bool calculate_observables       = false;
    bool check_causality            = false;
    fs::path results_directory        = "";

    static constexpr int VERBOSE      = 5;

    //==========================================================================
    // hydrodynamics settings

    // maximum upper limit for t
    double max_tau = 50.;

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

    bool perform_eos_checks           = false; ///< Flag to enable closure tests //TODO: Read its value from YAML file

    //==========================================================================
    // input parameter settings
    // quantities read in from InputParameters.inp file


    //------------------------------------
    // initial conditions
    string IC_type              = "";    ///< Type of initial conditions ("ccake", "iccing")
    string coordinate_system     = "hyperbolic"; ///< Coordinate system used in the initial conditions
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
    int rk_order                = 2;     /// Runge-Kutta order
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
    bool particlization_enabled = true;   ///< whether to use the "old" or "new"
                                          ///  particlization scheme
    double Freeze_Out_Temperature = 0.0;  ///< freeze-out temp. (at zero mu)

    //------------------------------------
    // equation of state
    string eos_type             = "";     ///< Type of equation of state ("conformal" or "table")
    string eos_path             = "";     ///< If "table", path to the equation of state file
    bool online_inverter_enabled = false;
    fs::path preinverted_eos_path;       ///< Path to the preinverted EOS file

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
    string etaMode              = "";    ///< Choose parametrization for eta/s
    double constant_eta_over_s  = 0.0;
    string shearRelaxMode       = "";    ///< Choose parametrization for relaxation time of shear mode

    //  - bulk quantities
    string zetaMode             = "";    ///< Choose parametrization for zeta/s
    double constant_zeta_over_s = 0.0;
    double cs2_dependent_zeta_A = 0.0;
    double cs2_dependent_zeta_p = 0.0;
    bool bulk_from_trace        = false; ///< Whether to use the trace of Tmunu to
                                          ///  calculate the bulk viscosity
    string bulkRelaxMode        = "";   ///< Choose parametrization for relaxation time of bulk mode
    bool modulate_zeta_with_tanh = true;   // forces zeta/s to decrease
                                           // smoothly to zero below
                                           // transition temperature
    // -- diffusion quantities
    string diffusionMode        = "";   ///< Choose parametrization for diffusion
    std::array<std::array<double, 3>, 3> kappa_matrix = {{{0.0, 0.0, 0.0},
                                                          {0.0, 0.0, 0.0},
                                                          {0.0, 0.0, 0.0}}};
    bool input_initial_diffusion = false; ///< Whether to use the initial diffusion
                                          ///  coefficients from the input file


    // make sure that all chosen settings make reasonable sense
    ///TODO: Check these are correct
    void check_consistency()
    {
      formatted_output::update("Impose consistency checks");

      // if eta/s == 0 identically, set using_shear to false
      if ( etaMode == "constant" && constant_eta_over_s < 1e-10 )
        using_shear  = false;
      else
        using_shear  = true;
      std::cout << "diffusionMode: " << diffusionMode << std::endl;
      if( diffusionMode != "disabled"){
        using_diffusion = true;
        std::cout << "using_diffusion: " << using_diffusion << std::endl;}
      else
        using_diffusion = false;
      return;
    }

};

typedef std::shared_ptr<Settings> pSettings;  // smart pointer to settings object

#endif