#ifndef INPUT_H
#define INPUT_H

#include <string>
#include <memory>
#include <filesystem>

#include <yaml-cpp/yaml.h>

#include "settings.h"
#include "defaults.h"
#include "formatted_output.h"

namespace fs = std::filesystem;
/// \brief Class for parsing input files and arguments

namespace ccake{
class Input
{
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