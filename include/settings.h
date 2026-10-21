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
    double shear_eta_0          = 0.08;  ///< eta/s at muB = 0   GeV (piecewise_muB mode)
    double shear_eta_2          = 0.08;  ///< eta/s at muB = 0.2 GeV (piecewise_muB mode)
    double shear_eta_4          = 0.08;  ///< eta/s at muB = 0.4 GeV (piecewise_muB mode)
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
    /// Use the EoS table's chiBB_nc field (a smooth non-critical baseline fitted to the same
    /// table) as the LOCAL xi_0 in the bulk critical scaling, instead of the global constant.
    /// Requires a 22-column table; ignored otherwise.  With this on the factor is normalised so
    /// that it is 1 away from the transition.
    bool critical_scaling_local_xi0 = false;
    /// Hard cap on the bulk critical enhancement factor. <= 0 disables the cap.  Needed because
    /// the raw factor reaches ~233 at the critical point and ~1259 at the plane maximum.
    double critical_scaling_max = 10.0;
    double constant_zeta_over_s = 0.0;
    double cs2_dependent_zeta_A = 0.0;
    double cs2_dependent_zeta_p = 0.0;
    double bulk_tau_pi_coeff    = 5.0; ///< prefactor in default bulk relaxation time tau_Pi = coeff*zeta/((1/3-cs2)^p*w)
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
    ///< path to the tabulated holographic baryon DC conductivity sigma_B(T,muB)
    string holographic_sigmaB_table_path = "";
    ///< dimensionless overall multiplier on the tabulated kappa_B = T*sigma_B
    double holographic_diffusion_scale = 1.0;
    ///< propagation speed squared of the diffusion front (nernst_einstein tau_q);
    ///< <= 0 means "use the local cs2"
    double diffusion_wave_speed_sq = 0.0;
    double b_q = 20.0;  ///< coefficient for tau_q = b_q * kappa_BB/n_B (kappa_over_n)
    bool tau_q_causality_floor = true; ///< enforce tau_q >= D_B (Kapusta-Young)
    ///< clamps on the nernst_einstein tau_q, in units of 1/T
    double tau_q_min_over_T = 0.1;
    double tau_q_max_over_T = 10.0;
    string diffusionRelaxMode = "default"; ///< Choose parametrization for diffusion relaxation time
    bool critical_scaling_diffusion = false; ///< Whether to enable Gaussian suppression of kappa near critical point
    double critical_point_T   = 0.0; ///< Critical point temperature [MeV]
    double critical_point_muB = 0.0; ///< Critical point baryon chemical potential [MeV]
    double critical_gaussian_width_T   = 0.0; ///< Gaussian width in T direction [MeV]
    double critical_gaussian_width_muB = 0.0; ///< Gaussian width in muB direction [MeV]
    /// Diffusion critical-scaling form. "gaussian_suppression" (legacy): kappa_B x [1 - Gaussian] centred on
    /// critical_point. "xi_power": the bulk scaling's form g(p) = min{[1 + (xi/xi_0)^p]/2, critical_scaling_max}
    string critical_scaling_diffusion_mode = "gaussian_suppression";
    double critical_scaling_kappa_power = 1.0;  ///< kappa_n ~ xi (model H: kappa_n ~ lambda_T ~ xi)
    double critical_scaling_tauq_power  = 2.0;  ///< tau_n ~ xi^2 (diffusive-shear slow mode with shear present)

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
    double source_test_particles = 1.0;   ///< divide source energy+charge deposits by this (SMASH N_test)
    bool baryon_source         = false;  ///< Whether to include baryon source term
    bool strangeness_source    = false;  ///< Whether to include strangeness source term
    bool electric_source       = false;  ///< Whether to include electric source term
    double smearing_radius     = 0.0;   ///< Transverse smearing radius for source terms [fm]
    double smearing_radius_eta = 0.0;   ///< Longitudinal (eta) smearing radius for source terms; used by the D=3 split deposit kernel
    string source_input_file   = "";    ///< Path to the source term file
    bool   source_propagate    = false; ///< If true, IC sources stream and lose energy at a constant rate; else instant dump at formation
    double source_loss_coefficient = 0.0; ///< propagate mode: constant loss rate, dE = loss_coefficient*dt [fm^-2]
    int    source_deposit_steps = 1;    ///< Instant mode: spread each dump over N steps (E0/N per step) to keep v<c; 1 = all at once
    bool   source_kernel_renormalize = false; ///< Shepard-renormalize each deposit by the discrete kernel sum (exact capture where fluid exists); OFF = legacy behavior
    double source_min_coverage = 0.1;   ///< kernel_renormalize only: defer a source's deposit while its kernel coverage sum is below this (fluid hasn't reached it yet)
    bool   source_reheat_frozen_at_birth = false; ///< particles below e_fo at tau0 are marked Freeze=5 (frozen at birth, no surface element yet) and return to Freeze=0 if a source deposit lifts them above e_fo; OFF = legacy (Freeze=4, permanent)
    bool   source_deposit_at_rapidity = false; ///< Deposit each source at eta_s = its momentum rapidity instead of its space-time rapidity, keeping the formation proper time (dynamical-initialisation prescription, e.g. 3D-Glauber+MUSIC places string sources at the string rapidity). At low energy the partons are made near z=0 while carrying |y| up to 2, and free-streaming cannot bridge that within the fireball lifetime
    double source_max_renormalization_factor = 2.0; ///< kernel_renormalize only: cap on the Shepard factor 1/coverage; below coverage=1/max_renormalization_factor only the safe fraction of the chunk lands and the budget carries forward

    //  - jet terms

    string jets_type            = "";    ///< Choose jet type and quantity
    int jets_Energy_scaling     = 0;    ///< Choose BBMG energy dependence
    int jets_Length_scaling     = 0;    ///< Choose BBMG length dependence
    int jets_Fluctuations       = 0;    ///< Choose BBMG fluctuations
    // 3D jet initialization mode:
    //   "file"      → read N jets from jets_input_file (cols: x y eta pT phi rapidity [pdg])
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