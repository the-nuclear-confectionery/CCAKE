#ifndef DEFAULTS_H
#define DEFAULTS_H

#include <map>
#include <string>
#include <vector>

typedef std::map  <std::string, std::string> setting_map;
typedef std::pair <std::string, std::string> setting_pair;

namespace ccake
{
  namespace defaults
  {
    const bool        online_inverter_enabled  = false;
    const bool        restrict_mu_T_ratios     = false;
    const bool        prohibit_unstable_cs2    = true;
    const bool        prohibit_acausal_cs2     = true;
    // gap-region fallback
    const std::string gap_table_path           = "";
    const std::string gap_lookup_mode          = "linear";
    const bool        gap_analytic_enabled     = false;
    const double      gap_match_tolerance      = 0.05;
    const double      max_tau                 = 50.;
    const bool        input_as_entropy        = false;
    const double      dt                      = 0.05;
    const double      hT                      = 0.3;
    const double      hEta                    = 0.3;
    const int         rk_order                = 2;
    const std::string kernel_type             = "cubic_spline";
    const double      e_cutoff                = 0.15;
    const bool        buffer_event            = false;
    const bool        circular_buffer         = false;
    const double      padding_thickness       = 0.1;
    const std::string eos_type                = "conformal";
    const std::string coordinate_system        = "hyperbolic";
    const bool        particlization_enabled  = true;
    const double      Freeze_Out_Temperature  = 150.; // MeV
    const std::string Freeze_Out_Type         = "fixed_T";
    const bool        baryon_charge_enabled   = true;
    const bool        strange_charge_enabled  = true;
    const bool        electric_charge_enabled = true;
    const std::string etaMode                 = "constant";
    const double      constant_eta_over_s     = 0.08;
    const std::string shearRelaxMode          = "default";
    const bool        input_initial_shear = false;
    const bool        use_vorticity        = false;
    const std::string delta_pipi_mode         = "israel-stewart";
    const std::string tau_pipi_mode          = "disabled";
    const std::string lambda_piPi_mode      = "disabled";
    const std::string phi6_mode             = "disabled";
    const std::string phi7_mode             = "disabled";
    const std::string zetaMode                = "cs2_dependent";
    const double      constant_zeta_over_s    = 0.005;
    const double      cs2_dependent_zeta_A    = 1.67552;  // 8*pi/15
    const double      cs2_dependent_zeta_p    = 2.0;
    const bool       critical_scaling_bulk   = false;
    const std::string bulkRelaxMode           = "default";
    const bool        modulate_zeta_with_tanh = true;
    const bool        delta_PiPi_mode       = "israel-stewart";
    const bool        lambda_Pipi_mode     = "disabled";
    const bool        phi1_mode          = "disabled";
    const bool        phi3_mode          = "disabled";
    const std::string diffusionMode           = "disabled";
    const std::string diffusionRelaxMode          = "default";
    const bool        critical_scaling_diffusion  = false;
    const double      critical_point_T            = 150.0;  // MeV
    const double      critical_point_muB          = 350.0;  // MeV
    const double      critical_gaussian_width_T   = 30.0;   // MeV
    const double      critical_gaussian_width_muB = 50.0;   // MeV
    const bool        input_initial_diffusion = false;
    const std::array<std::array<double, 3>, 3> 
                              constant_kappa_over_T2 = {0.0, 0.0, 0.0,
                                              0.0, 0.0, 0.0,
                                              0.0, 0.0, 0.0};
    const std::string relaxation_mode = "default";
    const std::string source_type            = "disabled";
    const std::string source_model           = "disabled";
    const double      source_normalization   = 1.0;
    const bool        baryon_source           = false;
    const bool        strangeness_source      = false;
    const bool        electric_source           = false;
    const double      smearing_radius        = 0.5;
    const double      smearing_radius_eta    = 0.5;
    const std::string source_input_file   = "disabled";
    const bool        source_propagate         = false;
    const double      source_loss_coefficient = 0.0;
    const double      source_entropy_ref       = 1.0;
    const int         source_deposit_steps     = 1;
    const std::string jets_type              = "disabled";
    const int         jets_Energy_scaling     = 0;
    const int         jets_Length_scaling     = 1;
    const int         jets_Fluctuations       = 0;
    const std::string jets_input_mode         = "hardcoded";
    const std::string jets_input_file         = "";
    const double      jets_pT                 = 10.0;
    const double      jets_phi                = 0.0;
    const double      jets_rapidity              = 0.0;
    const double      jets_x0                 = 0.0;
    const double      jets_y0                 = 0.0;
    const double      jets_eta0               = 0.0;
    const int         jets_njet               = 200000;
    const bool        jets_sample_rapidity       = false;
    const double      jets_rapidity_max          = 1.0;
    const bool        jets_auto_partner          = true;
    const bool        print_conservation_status = true;
    const bool        calculate_observables        = false;
    const bool        get_neighbors             = false;
    const bool        hdf_evolution            = false;
    const bool        txt_evolution            = false;
    const bool        jet_evolution            = false;
    const bool        check_causality         = false;
    const int         evolution_stride        = 1;
    const bool        causality_minimal       = false;
    const int         causality_minimal_stride = 20;
    const bool        bulk_from_trace       = false;
    // diffusion regulator (PhysRevC.98.034916, Eqs. C6-C8)
    const bool   diffusion_regulator_enabled  = false;
    const double diffusion_regulator_chi0     = 10.0;   ///< overall regulation strength [dimensionless]
    const double diffusion_regulator_e0       = 0.1;    ///< critical energy density e_0  [GeV/fm^3]
    const double diffusion_regulator_xi0      = 0.01;   ///< width parameter xi_0          [GeV/fm^3]
    const double diffusion_regulator_rq_max   = 1.0;    ///< maximum allowed r_q            [dimensionless]
  }
}

namespace parameter_settings
{
  inline setting_map get_defaults()
  {
    // set parameter defaults here
    std::vector<setting_pair>
      default_pairs
        = {
            setting_pair("ICtype",                 "ICCING"),
            setting_pair("ICfile",                 "initial_conditions/"
                                                   "Iccing_conditions.dat"),
            setting_pair("hT",                      "0.300000"),
            setting_pair("heta",                    "0.1"),
            setting_pair("rk_order",                "2"),
            setting_pair("dt",                     "0.05"),
            setting_pair("t0",                     "1.100000"),
            setting_pair("e_cutoff",               "0.00301"),
            setting_pair("EoS_Type",               "conformal"),
            setting_pair("EoS_Path",               "default"),
            setting_pair("etaMode",                "constant"),
            setting_pair("constant_eta_over_s",    "0.20"),
            setting_pair("shearRelaxMode",         "default"),
            setting_pair("zetaMode",               "constant"),
            setting_pair("constant_zeta_over_s",   "0.005"),
            setting_pair("cs2_dependent_zeta_A",   "1.67552"),  // 8*pi/15
            setting_pair("cs2_dependent_zeta_p",   "2.0"),
            setting_pair("bulkRelaxMode",          "default"),
            setting_pair("diffusionMode",          "disabled"),
            setting_pair("freezeoutT",             "150.000000"),
            setting_pair("freezeout",              "No_Freezeout"),
            setting_pair("Gubser_BSQmode",         "off")
          };

    // build defaults map and return it
    setting_map defaults;
    for ( auto & default_pair : default_pairs ) defaults.insert( default_pair );

    return defaults;
  }
}

#endif
