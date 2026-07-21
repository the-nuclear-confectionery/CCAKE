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
    bool using_vorticity             = false;  // whether to use vorticity in the shear
                                                // relaxation time calculation
    bool initializing_with_full_Tmunu = false;  // whether to initialize Pi from
                                                // tmunu_trace - p or not
    //==========================================================================
    // I/O settings
    bool txt_evolution              = false;
    bool hdf_evolution              = false;
    bool jet_evolution              = false;
    bool print_conservation_status  = true;
    bool get_neighbors              = false;
    bool calculate_observables       = false;
    bool check_causality            = false;
    int  evolution_stride           = 1;   ///< write full evolution output every Nth timestep (1 = every step)
    bool causality_minimal          = false; ///< dump per-particle causality diagnostics (x,y,eta,e,T,causality,invRe_sh/bk,Kn_sh/bk)
    int  causality_minimal_stride   = 20;  ///< minimal-causality dump every Nth timestep
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
    // freeze-in initial conditions (IC_type == "freezein"): read a raw
    // contravariant T^{mu nu} grid (AMPTGenesis raw_tmunu output) and match it
    // to ideal hydro fields with CCAKE's EoS (energy-momentum conserving,
    // arXiv:1103.4605 Eq. 14).
    fs::path freezein_input_file = "";    ///< Path to the raw T^{mu nu} grid file
    double freezein_tol          = 1e-10; ///< Convergence tol on the flow speed v
    int    freezein_max_iter     = 100;   ///< Max iterations of the matching loop
    bool   freezein_baryon_landau = false; ///< If true, n_a = j_a^mu u_mu (Landau scalar
                                           ///< rest density) instead of the flux-matched
                                           ///< j_a^tau/gamma. Keeps freeze-in energy/flow.

    //------------------------------------
    // simulations parameters
    double dt                   = 0.0;   ///< size of timestep
    double hT                   = 0.0;   ///< SPH kernel scale in transverse direction [fm]
    double hEta                 = 0.0;   ///< SPH kernel scale in longitudinal direction
    int rk_order                = 2;     /// Runge-Kutta order
    std::string kernel_type     = "";    ///< Which SPH kernel to use (cubic, quartic, quintic).
                                         ///  Currently only cubic is supported
                                         //TODO: Implement other kernels
    bool factorized_kernel      = false; ///< Use W_2D(r_perp,hT)*W_1D(deta,hEta) for D=3 SPH kernel
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
    bool normalize_by_T         = true;  ///< If true, table fields are dimensionless ratios (e.g. p/T^4);
    bool restrict_mu_T_ratios   = false; ///< If true, skip EoS solutions where max(|muB|,|muS|,|muQ|) > 3.5*T
                                         ///  if false, fields are in MeV^n and are rescaled by (1/hbarc)^n
    bool prohibit_unstable_cs2  = true;  ///< If true, reject EoS solutions with cs2 < 0 (spinodal/unstable branch).
                                         ///  Set false when using a table that intentionally includes the unstable phase.
    bool prohibit_acausal_cs2   = true;  ///< If true, reject EoS solutions with cs2 > 1 (acausal).

    // --- gap-region fallback ---
    fs::path    gap_table_path;              ///< optional scattered-points fallback file for the spinodal gap
    std::string gap_lookup_mode  = "linear"; ///< "nearest" or "linear" interpolation between gap points
    bool        gap_analytic_enabled = false;///< if true, try hard-coded analytic override before file lookup
    double      gap_match_tolerance  = 0.05; ///< max distance in (s,rhoB,rhoS,rhoQ) space to accept a gap hit

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
    bool use_vorticity          = false; ///< Whether to use vorticity in the shear
                                          ///  relaxation time calculation
    //  - delta pi/pi quantities
    string delta_pipi_mode      = "";    ///< Choose parametrization for delta pi/pi
    string tau_pipi_mode        = "";    ///< Choose parametrization for tau pi/pi
    string lambda_piPi_mode     = "";    ///< Choose parametrization for lambda pi/pi
    string phi6_mode            = "";    ///< Choose parametrization for phi6/pi  
    string phi7_mode            = "";    ///< Choose parametrization for phi7/pi
    bool input_initial_shear    = false;
    //  - bulk quantities
    string zetaMode             = "";    ///< Choose parametrization for zeta/s
    bool critical_scaling_bulk  = false; ///< Whether to enable critical scaling of bulk viscosity
    double constant_zeta_over_s = 0.0;
    double cs2_dependent_zeta_A = 0.0;
    double cs2_dependent_zeta_p = 0.0;
    bool bulk_from_trace        = false; ///< Whether to use the trace of Tmunu to
                                          ///  calculate the bulk viscosity
    string bulkRelaxMode        = "";   ///< Choose parametrization for relaxation time of bulk mode
    bool modulate_zeta_with_tanh = true;   // forces zeta/s to decrease
                                           // smoothly to zero below
                                           // transition temperature
    //  - other transport coefficients
    string delta_PiPi_mode    = "";   ///< Choose parametrization for delta Pi/Pi
    string lambda_Pipi_mode   = "";   ///< Choose parametrization for lambda Pi/Pi
    string phi1_mode          = "";   ///< Choose parametrization for phi1/Pi
    string phi3_mode          = "";   ///< Choose parametrization for phi3/Pi
    // -- diffusion quantities
    string diffusionMode        = "";   ///< Choose parametrization for diffusion
    std::array<std::array<double, 3>, 3> constant_kappa_over_T2 = {{{0.0, 0.0, 0.0},
                                                          {0.0, 0.0, 0.0},
                                                          {0.0, 0.0, 0.0}}};
    bool input_initial_diffusion = false; ///< Whether to use the initial diffusion
                                          ///  coefficients from the input file
    double C_B = 0.0;   ///< coefficient for coth-form baryon diffusion kappa_B (Eq. 16)
    string diffusionRelaxMode = "default"; ///< Choose parametrization for diffusion relaxation time
    bool critical_scaling_diffusion = false; ///< Whether to enable Gaussian suppression of kappa near critical point
    double critical_point_T   = 0.0; ///< Critical point temperature [MeV]
    double critical_point_muB = 0.0; ///< Critical point baryon chemical potential [MeV]
    double critical_gaussian_width_T   = 0.0; ///< Gaussian width in T direction [MeV]
    double critical_gaussian_width_muB = 0.0; ///< Gaussian width in muB direction [MeV]

    // -- diffusion regulator (PhysRevC.98.034916, Appendix C, Eqs. C6-C8)
    bool   diffusion_regulator_enabled  = false; ///< enable/disable the diffusion regulator
    double diffusion_regulator_chi0     = 10.0;  ///< overall regulation strength chi_0 [dimensionless]
    double diffusion_regulator_e0       = 0.1;   ///< critical energy density e_0  [GeV/fm^3]
    double diffusion_regulator_xi0      = 0.01;  ///< width parameter xi_0          [GeV/fm^3]
    double diffusion_regulator_rq_max   = 1.0;   ///< maximum allowed r_q value      [dimensionless]

    //  - source terms

    string source_type           = "";    ///< Choose source term type for shear
    string source_model          = "";    ///< Choose source term model
    double source_normalization  = 1.0;   ///< Normalization factor for source terms
    bool baryon_source         = false;  ///< Whether to include baryon source term
    bool strangeness_source    = false;  ///< Whether to include strangeness source term
    bool electric_source       = false;  ///< Whether to include electric source term
    double smearing_radius     = 0.0;   ///< Transverse smearing radius for source terms [fm]
    double smearing_radius_eta = 0.0;   ///< Longitudinal (eta) smearing radius for source terms; used by the D=3 split deposit kernel
    string source_input_file   = "";    ///< Path to the source term file
    bool   source_propagate    = false; ///< If true, SMASH sources propagate & deplete like a jet; else instant dump at t_src
    double source_loss_coefficient = 0.0; ///< CTE in dE = CTE * s(x)/s_0 * dt for propagating SMASH sources
    double source_entropy_ref      = 1.0; ///< s_0 reference entropy density for propagating SMASH sources
    int    source_deposit_steps = 1;    ///< Instant mode: spread each dump over N steps (E0/N per step) to keep v<c; 1 = all at once

    //  - jet terms

    string jets_type            = "";    ///< Choose jet type and quantity
    int jets_Energy_scaling     = 0;    ///< Choose BBMG energy dependence
    int jets_Length_scaling     = 0;    ///< Choose BBMG length dependence
    int jets_Fluctuations       = 0;    ///< Choose BBMG fluctuations
    // 3D jet initialization mode:
    //   "file"      → read N jets from jets_input_file (cols: x y eta pT phi rapidity)
    //   "hardcoded" → single jet at (x0,y0,eta0) with (pT, phi, rapidity)
    //   "sampled"   → sample njet jets: positions from SPH particles above Freezeout_Temp,
    //                 phi from phi_bins or uniform, rapidity sampled per sampling block.
    string jets_input_mode      = "hardcoded";
    string jets_input_file      = "";
    double jets_pT              = 0.0;   ///< Jet transverse momentum [GeV]
    double jets_phi             = 0.0;   ///< Jet azimuthal angle [rad]
    double jets_rapidity           = 0.0;   ///< Jet rapidity
    double jets_x0              = 0.0;   ///< Jet initial x [fm]
    double jets_y0              = 0.0;   ///< Jet initial y [fm]
    double jets_eta0            = 0.0;   ///< Jet initial spacetime rapidity
    // Sampling subsection (only used when input_mode == "sampled").
    // phi is always sampled uniformly in [0, 2π) in this mode.
    int    jets_njet            = 200000; ///< Number of jets per event to sample
    bool   jets_sample_rapidity    = false;  ///< true → uniform in [-rapidity_max, rapidity_max]; false → fixed at jets_rapidity
    double jets_rapidity_max       = 1.0;    ///< Half-width of jet rapidity sampling range
    // file mode only: true → auto back-to-back partner (φ+π, same rapidity) per line;
    //                 false → one jet per line (faithful explicit dijet/parton list).
    bool   jets_auto_partner    = true;

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
      if ( use_vorticity == true){
        using_vorticity = true;
      }
      else
        using_vorticity = false;

      return;
    }

};

typedef std::shared_ptr<Settings> pSettings;  // smart pointer to settings object

#endif