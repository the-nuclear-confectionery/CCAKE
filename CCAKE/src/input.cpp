#ifndef INPUT_CPP
#define INPUT_CPP

#include "input.h"

using namespace ccake;

/// \brief Constructor of the Input class. Parses the input file and sets the
/// values of the settings object.
/// \param[in] argc Number of command line arguments
/// \param[in] argv Array of command line arguments
Input::Input( int argc, char** argv)
{
  check_args( argc, argv ); //Will abort execution if arguments are invalid
  load_settings_file();
  ////read_in_initial_conditions();
  ////validate_input();


  formatted_output::report( "Input parameters file: "
                            + settings_file_path.string() );
  formatted_output::report( "All results will be stored in: "
                            + results_directory.string() );
}

/// \brief Check if input arguments are valid
/// \param[in] argc Number of command line arguments
/// \param[in] argv Array of command line arguments
void Input::check_args( int argc, char** argv )
{
  //Check if right number of arguments were fed in. If not, print usage
  ///and abort execution.
  if ( argc != 3 )
  {
    formatted_output::report("Error: Wrong number of arguments!");
    std::cerr << "Usage: " << argv[0] << " <path_to_settings_file> <path_to_results_directory>" << endl;
    std::cerr << "       Please cite all our papers." << std::endl;
    exit(EXIT_FAILURE);
  }

  //Check if settings file exists. If not, print error message and abort execution.
  settings_file_path = argv[1];
  if ( !fs::exists(settings_file_path) )
  {
      std::cout << "Error: " << settings_file_path << " does not exist!" << endl;
      exit(EXIT_FAILURE);
  }

  results_directory = argv[2];
  if ( !fs::exists(results_directory) )
  {
    fs::create_directory(results_directory);
  }

}


/// \brief Reads in the settings file. Make sure to initialize the pointer to
/// the settings object before calling this function.
/// \param[in] path_to_settings_file Path to the settings file.
void Input::load_settings_file()
{
  formatted_output::report("Reading in input parameter settings");
  settingsPtr = std::make_shared<Settings>();
  YAML::Node input_file = YAML::LoadFile(settings_file_path.string());
  if ( !decode_settings(input_file) )
  {
    formatted_output::report("Error: Could not decode settings file!");
    exit(EXIT_FAILURE);
  }
  return;
}

bool Input::decode_settings(const YAML::Node& node){

    //Initial conditions node
    settingsPtr->IC_type = node["initial_conditions"]["type"].as<std::string>();
    settingsPtr->IC_file = fs::path(node["initial_conditions"]["file"].as<std::string>());
    settingsPtr->t0 = node["initial_conditions"]["t0"].as<double>();
    settingsPtr->dim = node["initial_conditions"]["dimension"].as<unsigned int>();
    settingsPtr->input_as_entropy = node["initial_conditions"]["input_as_entropy"].as<bool>();

    //Parameters node
    settingsPtr->dt = node["parameters"]["dt"].as<double>();
    settingsPtr->hT = node["parameters"]["h_T"].as<double>();
    settingsPtr->hEta = node["parameters"]["h_eta"].as<double>();
    settingsPtr->kernel_type = node["parameters"]["kernel_type"].as<std::string>();
    settingsPtr->e_cutoff = node["parameters"]["energy_cutoff"].as<double>();

    //buffer_particles subnode
    settingsPtr->buffer_event = node["parameters"]["buffer_particles"]["enabled"].as<bool>();
    settingsPtr->circular_buffer = node["parameters"]["buffer_particles"]["circular"].as<bool>();
    settingsPtr->padding_thickness = node["parameters"]["buffer_particles"]["padding_thickness"].as<double>();

    //eos node
    settingsPtr->eos_type = node["eos"]["type"].as<std::string>();
    settingsPtr->eos_path = fs::path(node["eos"]["path"].as<std::string>());

    //particlization node
    settingsPtr->particlization_enabled = node["particlization"]["enabled"].as<bool>();
    settingsPtr->Freeze_Out_Temperature = node["particlization"]["T"].as<double>();
    settingsPtr->Freeze_Out_Type = node["particlization"]["type"].as<string>();

    //hydro node
    settingsPtr->baryon_charge_enabled = node["hydro"]["baryon_charge_enabled"].as<bool>();
    settingsPtr->strange_charge_enabled = node["hydro"]["strange_charge_enabled"].as<bool>();
    settingsPtr->electric_charge_enabled = node["hydro"]["electric_charge_enabled"].as<bool>();

    //viscous_parameters subnode
    //shear subsubnode
    settingsPtr->etaMode = node["hydro"]["viscous_parameters"]["shear"]["mode"].as<std::string>();
    settingsPtr->constant_eta_over_s = node["hydro"]["viscous_parameters"]["shear"]["constant_eta_over_s"].as<double>();
    settingsPtr->shearRelaxMode = node["hydro"]["viscous_parameters"]["shear"]["relaxation_mode"].as<std::string>();
    //bulk subsubnode
    settingsPtr->zetaMode = node["hydro"]["viscous_parameters"]["bulk"]["mode"].as<std::string>();
    settingsPtr->constant_zeta_over_s = node["hydro"]["viscous_parameters"]["bulk"]["constant_zeta_over_s"].as<double>();
    settingsPtr->bulkRelaxMode = node["hydro"]["viscous_parameters"]["bulk"]["relaxation_mode"].as<std::string>();
    settingsPtr->cs2_dependent_zeta_A = node["hydro"]["viscous_parameters"]["bulk"]["cs2_dependent_zeta_A"].as<double>();
    settingsPtr->cs2_dependent_zeta_p = node["hydro"]["viscous_parameters"]["bulk"]["cs2_dependent_zeta_p"].as<double>();
    settingsPtr->modulate_zeta_with_tanh = node["hydro"]["viscous_parameters"]["bulk"]["modulate_with_tanh"].as<bool>();

    return true;

}

#endif