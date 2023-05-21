#ifndef INPUT_CPP
#define INPUT_CPP

#include "input.h"

namespace cc = ccake;
namespace consts = constants;
//using namespace ccake;

/// \brief Constructor of the Input class. Parses the input file and sets the
/// values of the settings object.
/// \param[in] argc Number of command line arguments
/// \param[in] argv Array of command line arguments
cc::Input::Input( int argc, char** argv)
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

/// \brief Check if input arguments are valid. If not, print usage and abort.
/// This creates the output directory if it does not exist.
/// \param[in] argc Number of command line arguments
/// \param[in] argv Array of command line arguments
void cc::Input::check_args( int argc, char** argv )
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


/// \brief Parses the settings file. You need to have checked the input
/// arguments first.
/// \param[in] path_to_settings_file Path to the settings file.
void cc::Input::load_settings_file()
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

/// \brief Decodes the settings file.
/// \param[in] node YAML root node containing the settings file.
/// \return True if successful, false otherwise.
bool cc::Input::decode_settings(const YAML::Node& node){
    //TODO: Maybe its possible to coalesce this into a loop?
    //--------------------------------------------------------------------------
    //Initial conditions node
    try{
      settingsPtr->IC_type = node["initial_conditions"]["type"].as<std::string>();
    } catch (...){
      formatted_output::detail("ERROR: Could not read initial conditions type!");
      formatted_output::detail("This is a mandatory parameter. Aborting execution.");
      return false;
    }

    try{
      settingsPtr->IC_file = fs::path(node["initial_conditions"]["file"].as<std::string>());
    } catch (...) {
      formatted_output::detail("ERROR: Could not read initial conditions file!");
      formatted_output::detail("This is a mandatory parameter. Aborting execution.");
      return false;
    }

    try {
      settingsPtr->t0 = node["initial_conditions"]["t0"].as<double>();
    }
    catch (...){
      formatted_output::detail("ERROR: Could not read initial conditions intiyal time t0!");
      formatted_output::detail("This is a mandatory parameter. Aborting execution.");
      return false;
    }

    try {
      settingsPtr->dim = node["initial_conditions"]["dimension"].as<unsigned int>();
    } catch (...) {
      formatted_output::detail("ERROR: Could not read initial conditions dimension!");
      formatted_output::detail("This is a mandatory parameter. Aborting execution.");
      return false;
    }

    try {
      settingsPtr->input_as_entropy = node["initial_conditions"]["input_as_entropy"].as<bool>();
    } catch (...) {
      formatted_output::report("WARNING: Could not read initial condition input_as_entropy!");
      formatted_output::report("This is an optional parameter. Setting to false.");
      settingsPtr->input_as_entropy = cc::defaults::input_as_entropy;
    }

    //--------------------------------------------------------------------------
    //Parameters node
    try {
      settingsPtr->dt = node["parameters"]["dt"].as<double>();
    } catch (...) {
      formatted_output::report("WARNING: Could not read parameters dt!");
      formatted_output::report("This is an optional parameter. Setting to default value.");
      settingsPtr->dt = cc::defaults::dt;
    }

    try {
      settingsPtr->hT = node["parameters"]["h_T"].as<double>();
    } catch (...) {
      formatted_output::report("WARNING: Could not read parameters h_T!");
      formatted_output::report("This is an optional parameter. Setting to default value.");
      settingsPtr->hT = cc::defaults::hT;
    }
    
    try {
      settingsPtr->hEta = node["parameters"]["h_eta"].as<double>();
    } catch (...) {
      formatted_output::report("WARNING: Could not read parameters h_eta!");
      formatted_output::report("This is an optional parameter. Setting to default value.");
      settingsPtr->hEta = cc::defaults::hEta;
    }

    try {
      settingsPtr->kernel_type = node["parameters"]["kernel_type"].as<std::string>();
    } catch (...){
      formatted_output::report("WARNING: Could not read parameters kernel_type!");
      formatted_output::report("This is an optional parameter. Setting to default value.");
      settingsPtr->kernel_type = cc::defaults::kernel_type;
    }

    try {
      settingsPtr->e_cutoff = node["parameters"]["energy_cutoff"].as<double>();
    } catch (...) {
      formatted_output::report("WARNING: Could not read parameters energy_cutoff!");
      formatted_output::report("This is an optional parameter. Setting to default value.");
      settingsPtr->e_cutoff = cc::defaults::e_cutoff;
    }

    //--------------------------------------------------------------------------
    //buffer_particles subnode
    try {
      settingsPtr->buffer_event = node["parameters"]["buffer_particles"]["enabled"].as<bool>();
    } catch (...) {
      formatted_output::report("WARNING: Could not read parameters buffer_particles enabled!");
      formatted_output::report("This is an optional parameter. Setting to default value.");
      settingsPtr->buffer_event = cc::defaults::buffer_event;
    }

    if (settingsPtr->buffer_event) {
      try {
        settingsPtr->circular_buffer = node["parameters"]["buffer_particles"]["circular"].as<bool>();
      } catch (...) {
        formatted_output::report("WARNING: Could not read parameters buffer_particles circular!");
        formatted_output::report("This is an optional parameter. Setting to default value.");
        settingsPtr->circular_buffer = cc::defaults::circular_buffer;
      }

      try{
        settingsPtr->padding_thickness = node["parameters"]["buffer_particles"]["padding_thickness"].as<double>();
      } catch (...) {
        formatted_output::report("WARNING: Could not read parameters buffer_particles padding_thickness!");
        formatted_output::report("This is an optional parameter. Setting to default value.");
        settingsPtr->padding_thickness = cc::defaults::padding_thickness;
      }
    }

    //--------------------------------------------------------------------------
    //eos node
    try {
      settingsPtr->eos_type = node["eos"]["type"].as<std::string>();
    } catch (...) {
      formatted_output::detail("ERROR: Could not read eos type!");
      formatted_output::detail("This is an optional parameter. Setting to default value.");
      settingsPtr->eos_type = cc::defaults::eos_type;
    }

    if (settingsPtr->eos_type == "table"){
      try {
        settingsPtr->eos_path = fs::path(node["eos"]["path"].as<std::string>());
      } catch (...) {
        formatted_output::report("ERROR: Could not read eos path!");
        formatted_output::report("This is a mandatory parameter for table eos. Aborting execution,");
        return false;
      }
    }

    //--------------------------------------------------------------------------
    //particlization node
    try {
      settingsPtr->particlization_enabled = node["particlization"]["enabled"].as<bool>();
    } catch (...) {
      formatted_output::report("WARNING: Could not read particlization enabled!");
      formatted_output::report("This is an optional parameter. Setting to default value.");
      settingsPtr->particlization_enabled = cc::defaults::particlization_enabled;
    }

    if (settingsPtr->particlization_enabled){
      try {
        settingsPtr->Freeze_Out_Temperature = node["particlization"]["T"].as<double>();
      } catch (...) {
        formatted_output::report("WARNING: Could not read particlization T!");
        formatted_output::report("This is an optional parameter. Setting to default value.");
        settingsPtr->Freeze_Out_Temperature = cc::defaults::Freeze_Out_Temperature;
      }

      try {
        settingsPtr->Freeze_Out_Type = node["particlization"]["type"].as<string>();
      } catch (...) {
        formatted_output::report("WARNING: Could not read particlization type!");
        formatted_output::report("This is an optional parameter. Setting to default value.");
        settingsPtr->Freeze_Out_Type = cc::defaults::Freeze_Out_Type;
      }
    }

    //--------------------------------------------------------------------------
    //hydro node
    try {
      settingsPtr->baryon_charge_enabled = node["hydro"]["baryon_charge_enabled"].as<bool>();
    } catch (...) {
      formatted_output::report("WARNING: Could not read hydro baryon_charge_enabled!");
      formatted_output::report("This is an optional parameter. Setting to default value.");
      settingsPtr->baryon_charge_enabled = cc::defaults::baryon_charge_enabled;
    }

    try {
      settingsPtr->strange_charge_enabled = node["hydro"]["strange_charge_enabled"].as<bool>();
    } catch (...) {
      formatted_output::report("WARNING: Could not read hydro strange_charge_enabled!");
      formatted_output::report("This is an optional parameter. Setting to default value.");
      settingsPtr->strange_charge_enabled = cc::defaults::strange_charge_enabled;
    }

    try {
      settingsPtr->electric_charge_enabled = node["hydro"]["electric_charge_enabled"].as<bool>();
    } catch (...) {
      formatted_output::report("WARNING: Could not read hydro electric_charge_enabled!");
      formatted_output::report("This is an optional parameter. Setting to default value.");
      settingsPtr->electric_charge_enabled = cc::defaults::electric_charge_enabled;
    }
    //--------------------------------------------------------------------------
    //viscous_parameters subnode
    //shear subsubnode
    try {
      settingsPtr->etaMode = node["hydro"]["viscous_parameters"]["shear"]["mode"].as<std::string>();
    } catch (...) {
      formatted_output::report("WARNING: Could not read viscous_parameters shear mode!");
      formatted_output::report("This is an optional parameter. Setting to default value.");
      settingsPtr->etaMode = cc::defaults::etaMode;
    }

    if (settingsPtr->etaMode == "constant"){
      try{
        settingsPtr->constant_eta_over_s = node["hydro"]["viscous_parameters"]["shear"]["constant_eta_over_s"].as<double>();
      } catch (...) {
        formatted_output::report("WARNING: Could not read viscous_parameters shear constant_eta_over_s!");
        formatted_output::report("This is an optional parameter. Setting to default value.");
        settingsPtr->constant_eta_over_s = cc::defaults::constant_eta_over_s;
      }
    }

    try {
      settingsPtr->shearRelaxMode = node["hydro"]["viscous_parameters"]["shear"]["relaxation_mode"].as<std::string>();
    } catch (...) {
      formatted_output::report("WARNING: Could not read viscous_parameters shear relaxation_mode!");
      formatted_output::report("This is an optional parameter. Setting to default value.");
      settingsPtr->shearRelaxMode = cc::defaults::shearRelaxMode;
    }

    //bulk subsubnode
    try {
      settingsPtr->zetaMode = node["hydro"]["viscous_parameters"]["bulk"]["mode"].as<std::string>();
    } catch (...) {
      formatted_output::report("WARNING: Could not read viscous_parameters bulk mode!");
      formatted_output::report("This is an optional parameter. Setting to default value.");
      settingsPtr->zetaMode = cc::defaults::zetaMode;
    }

    if (settingsPtr->zetaMode == "constant"){
      try{
        settingsPtr->constant_zeta_over_s = node["hydro"]["viscous_parameters"]["bulk"]["constant_zeta_over_s"].as<double>();
      } catch (...) {
        formatted_output::report("WARNING: Could not read viscous_parameters bulk constant_zeta_over_s!");
        formatted_output::report("This is an optional parameter. Setting to default value.");
        settingsPtr->constant_zeta_over_s = cc::defaults::constant_zeta_over_s;
      }
    }
    if (settingsPtr->zetaMode == "cs2_dependent"){
      try {
        settingsPtr->cs2_dependent_zeta_A = node["hydro"]["viscous_parameters"]["bulk"]["cs2_dependent_zeta_A"].as<double>();
      } catch (...) {
        formatted_output::report("WARNING: Could not read viscous_parameters bulk cs2_dependent_zeta_A!");
        formatted_output::report("This is an optional parameter. Setting to default value.");
        settingsPtr->cs2_dependent_zeta_A = cc::defaults::cs2_dependent_zeta_A;
      }

      try {
        settingsPtr->cs2_dependent_zeta_p = node["hydro"]["viscous_parameters"]["bulk"]["cs2_dependent_zeta_p"].as<double>();
      } catch (...) {
        formatted_output::report("WARNING: Could not read viscous_parameters bulk cs2_dependent_zeta_p!");
        formatted_output::report("This is an optional parameter. Setting to default value.");
        settingsPtr->cs2_dependent_zeta_p = cc::defaults::cs2_dependent_zeta_p;
      }
    }

    try{
      settingsPtr->bulkRelaxMode = node["hydro"]["viscous_parameters"]["bulk"]["relaxation_mode"].as<std::string>();
    } catch (...) {
      formatted_output::report("WARNING: Could not read viscous_parameters bulk relaxation_mode!");
      formatted_output::report("This is an optional parameter. Setting to default value.");
      settingsPtr->bulkRelaxMode = cc::defaults::bulkRelaxMode;
    }

    try{
      settingsPtr->modulate_zeta_with_tanh = node["hydro"]["viscous_parameters"]["bulk"]["modulate_with_tanh"].as<bool>();
    } catch (...) {
      formatted_output::report("WARNING: Could not read viscous_parameters bulk modulate_with_tanh!");
      formatted_output::report("This is an optional parameter. Setting to default value.");
      settingsPtr->modulate_zeta_with_tanh = cc::defaults::modulate_zeta_with_tanh;
    }
    return true;
}

#endif