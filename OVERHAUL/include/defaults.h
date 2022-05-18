#ifndef DEFAULTS_H
#define DEFAULTS_H

#include <map>
#include <string>

typedef std::map  <std::string, std::string> setting_map;
typedef std::pair <std::string, std::string> setting_pair;


namespace parameter_settings
{
  setting_map get_defaults()
  {
    setting_map defaults;

    defaults.insert( setting_pair("ICtype",     "ICCING") );
    defaults.insert( setting_pair("ICoption",   "Default") );
    defaults.insert( setting_pair("ICfile",     "initial_conditions/"
                                                "Iccing_conditions.dat") );
    defaults.insert( setting_pair("h",          "0.300000") );
    defaults.insert( setting_pair("dt",         "0.05") );
    defaults.insert( setting_pair("t0",         "1.100000") );
    defaults.insert( setting_pair("EoS_Type",   "Conformal") );
    defaults.insert( setting_pair("EoS_Option", "Default") );
    defaults.insert( setting_pair("eta",        "constant") );
    defaults.insert( setting_pair("etaOpt",     "0.20") );
    defaults.insert( setting_pair("shearRelax", "Default") );
    defaults.insert( setting_pair("zeta",       "constant") );
    defaults.insert( setting_pair("zetaOpt",    "0.005") );
    defaults.insert( setting_pair("bulkRelax",  "Default") );
    defaults.insert( setting_pair("freezeoutT", "150.000000") );
    defaults.insert( setting_pair("freezeout",  "No_Freezeout") );

    return defaults;
  }
}

#endif