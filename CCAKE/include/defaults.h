#ifndef DEFAULTS_H
#define DEFAULTS_H

#include <map>
#include <string>
#include <vector>

typedef std::map  <std::string, std::string> setting_map;
typedef std::pair <std::string, std::string> setting_pair;


namespace parameter_settings
{
  inline setting_map get_defaults()
  {
    // set parameter defaults here
    std::vector<setting_pair>
      default_pairs
        = {
            setting_pair("ICtype",                 "ICCING"),
            setting_pair("ICoption",               "default"),
            setting_pair("ICfile",                 "initial_conditions/"
                                                   "Iccing_conditions.dat"),
            setting_pair("h",                      "0.300000"),
            setting_pair("dt",                     "0.05"),
            setting_pair("t0",                     "1.100000"),
            setting_pair("e_cutoff",               "0.00301"),
            setting_pair("EoS_Type",               "conformal"),
            setting_pair("EoS_Option",             "default"),
            setting_pair("etaMode",                "constant"),
            setting_pair("constant_eta_over_s",    "0.20"),
            setting_pair("shearRelaxMode",         "default"),
            setting_pair("zetaMode",               "constant"),
            setting_pair("constant_zeta_over_s",   "0.005"),
            setting_pair("cs2_dependent_zeta_A",   "1.67552"),  // 8*pi/15
            setting_pair("cs2_dependent_zeta_p",   "2.0"),
            setting_pair("bulkRelaxMode",          "default"),
            setting_pair("freezeoutT",             "150.000000"),
            setting_pair("freezeout",              "No_Freezeout")
          };
      
    // build defaults map and return it
    setting_map defaults;
    for ( auto & default_pair : default_pairs) defaults.insert( default_pair );

    return defaults;
  }
}

#endif