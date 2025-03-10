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
    const std::string zetaMode                = "cs2_dependent";
    const double      constant_zeta_over_s    = 0.005;
    const double      cs2_dependent_zeta_A    = 1.67552;  // 8*pi/15
    const double      cs2_dependent_zeta_p    = 2.0;
    const std::string bulkRelaxMode           = "default";
    const bool        modulate_zeta_with_tanh = true;
    const bool        print_conservation_status = true;
    const bool        calculate_observables        = false;
    const bool        hdf_evolution            = false;
    const bool        txt_evolution            = false;
    const bool        check_causality         = false;
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
