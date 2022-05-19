#ifndef DEFAULTS_BLAH_H
#define DEFAULTS_BLAH_H

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
            setting_pair("ICoption",               "Default"),
            setting_pair("ICfile",                 "initial_conditions/"
                                                   "Iccing_conditions.dat"),
            setting_pair("h",                      "0.300000"),
            setting_pair("dt",                     "0.05"),
            setting_pair("t0",                     "1.100000"),
            setting_pair("EoS_Type",               "Conformal"),
            setting_pair("EoS_Option",             "Default"),
            setting_pair("eta",                    "constant"),
            setting_pair("constant_eta_over_s",    "0.20"),
            setting_pair("shearRelax",             "Default"),
            setting_pair("zeta",                   "constant"),
            setting_pair("constant_zeta_over_s",   "0.005"),
            setting_pair("cs2_dependent_zeta_A",   "1.67552"),  // 8*pi/15
            setting_pair("cs2_dependent_zeta_p",   "2.0"),
            setting_pair("bulkRelax",              "Default"),
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