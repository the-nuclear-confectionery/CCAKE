#ifndef INPUT_H
#define INPUT_H

#include <string>
#include <memory>
#include <filesystem>

#include <yaml-cpp/yaml.h>
#include <Kokkos_Core.hpp>

#include "settings.h"
#include "defaults.h"
#include "formatted_output.h"

namespace fs = std::filesystem;
/// \brief Class for parsing input files and arguments

namespace ccake{
class Input
{
/*******************************************************************************
* The Input class is responsible for parsing the input file and command-line
* arguments. The constructor receives the command-line arguments, just as the
* main function does. `argv[1]` is the path to the input file, and `argv[2]` is
* the path to the results directory.
*
* We verify if the input file exists, and if it does, we parse it, initiallizing
* the settings object. The settings object is created in the heap and stored as
* a `shared pointer`. You should get ownership of the pointer before destroying
* the Input object. The settings object can then be shared with the other
* classes that needs it.
*
* We also create the results directory if it does not exist.
 ******************************************************************************/

private:
  void check_args( int argc, char** argv );
  bool decode_settings(const YAML::Node& node);
  void load_settings_file();
  void set_output_directory(std::string path_to_results_directory);
  //void read_in_initial_conditions();
  //void validate_input();


  //Helper functions for parsing input files
  //void set_value( setting_map & values, const std::string key,
  //                                      const std::string value );
  //std::string get_value( setting_map & values, const std::string & key );

  //Data-members
  fs::path settings_file_path; ///< Path to the settings file
  fs::path results_directory;  ///< Path to the results directory

public:
  Input() = delete; ///< Default constructor is deleted. You must specify a path
                    /// to the input file and the results directory.
  Input( int argc, char** argv );
  ~Input(){};
  std::shared_ptr<Settings> settingsPtr; ///< Pointer to the settings object

  //Getters
};
} // namespace ccake
#endif