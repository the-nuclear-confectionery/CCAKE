#include <filesystem>
#include <ostream>
#include <gtest/gtest.h>

#include "input.h"

//TODO: Write death tests for the case of missing fields in the input file

namespace fs = std::filesystem;
#define TOL 1.E-14

//Forward declarations of helper functions
void create_test_input();

/// \brief Happy scenartio: Create a sample input and check if the constructor
/// reads it correctly
TEST (InputTest, Constructor){
  create_test_input();
  fs::create_directory("test_result_output");

  std::array<std::string,3> fake_args = {"InputTest.Constructor", "sample_input.yaml", "test_result_output"};
  char** man_args;
  man_args = new char*[3];
  for(int i=0; i < 3; i++) {
    man_args[i] = fake_args[i].data();
  }

  ccake::Input in(fake_args.size(), man_args);

  //Check if the input file was read correctly

  //IC parameters
  EXPECT_EQ( in.settingsPtr->IC_type, "ccake" );
  EXPECT_EQ( in.settingsPtr->IC_file, "initial_conditions/ccake_trento.dat" );
  EXPECT_DOUBLE_EQ( in.settingsPtr->t0, 0.6);
  EXPECT_EQ( in.settingsPtr->dim, 2 );
  EXPECT_TRUE( in.settingsPtr->input_as_entropy );

  //Simulation parameters
  EXPECT_DOUBLE_EQ( in.settingsPtr->dt, 0.03 );
  EXPECT_DOUBLE_EQ( in.settingsPtr->hT, 0.5 );
  EXPECT_DOUBLE_EQ( in.settingsPtr->hEta, 0.5 );
  EXPECT_EQ( in.settingsPtr->kernel_type, "quartic_spline" );
  EXPECT_DOUBLE_EQ( in.settingsPtr->e_cutoff, .1 );
  EXPECT_TRUE( in.settingsPtr->buffer_event );
  EXPECT_TRUE (in.settingsPtr->circular_buffer);
  EXPECT_DOUBLE_EQ( in.settingsPtr->padding_thickness, .2);

  //EOS parameters
  EXPECT_EQ( in.settingsPtr->eos_type, "table" );
  EXPECT_EQ( in.settingsPtr->eos_path, "Houston" );

  //Particlization settings
  EXPECT_TRUE( in.settingsPtr->particlization_enabled);
  EXPECT_EQ( in.settingsPtr->Freeze_Out_Type, "fixed_T" );
  EXPECT_DOUBLE_EQ( in.settingsPtr->Freeze_Out_Temperature, 140. );

  //Hydro values
  EXPECT_FALSE( in.settingsPtr->baryon_charge_enabled );
  EXPECT_FALSE( in.settingsPtr->strange_charge_enabled );
  EXPECT_FALSE( in.settingsPtr->electric_charge_enabled );

  //Shear parameters
  EXPECT_EQ( in.settingsPtr->etaMode, "constant" );
  EXPECT_DOUBLE_EQ( in.settingsPtr->constant_eta_over_s, 0.2 );
  EXPECT_EQ( in.settingsPtr->shearRelaxMode, "default" );

  //Bulk parameters
  EXPECT_EQ( in.settingsPtr->zetaMode, "cs2_dependent" );
  EXPECT_DOUBLE_EQ( in.settingsPtr->constant_zeta_over_s, 0.0 );
  EXPECT_EQ( in.settingsPtr->bulkRelaxMode, "default" );
  EXPECT_DOUBLE_EQ( in.settingsPtr->cs2_dependent_zeta_A, 1.67552 );
  EXPECT_DOUBLE_EQ( in.settingsPtr->cs2_dependent_zeta_p, 2.0 );
  EXPECT_TRUE( in.settingsPtr->modulate_zeta_with_tanh);







  fs::remove_all("test_result_output");
  //fs::remove("sample_input.yaml");

}

/// \brief This is a scenario where the output directory does not exists. The
/// constructor should create it.
TEST (InputTest, ConstructorMkDir){
  create_test_input();

  std::array<std::string,3> fake_args = {"InputTest.Constructor", "sample_input.yaml", "test_result_output"};
  char** man_args;
  man_args = new char*[3];
  for(int i=0; i < 3; i++) {
    man_args[i] = fake_args[i].data();
  }

  ccake::Input in(fake_args.size(), man_args);


  EXPECT_TRUE( fs::exists("test_result_output"));

  fs::remove_all("test_result_output");
  //fs::remove("sample_input.yaml");

}

/// \brief This is a scenario where the input file does not exists. The program
/// should abort.
TEST (InputDeathTest, Constructor){

  fs::create_directory("test_result_output");

  std::array<std::string,3> fake_args = {"InputTest.Constructor", "invalid.yaml", "test_result_output"};
  char** man_args;
  man_args = new char*[3];
  for(int i=0; i < 3; i++) {
    man_args[i] = fake_args[i].data();
  }
  EXPECT_EXIT( ccake::Input(fake_args.size(), man_args) , testing::ExitedWithCode(EXIT_FAILURE), ".*");

}

/// \brief Helper function that generates a sample input file, as it should be
/// expected by the code
void create_test_input(){
  std::ofstream sample_input("sample_input.yaml",std::ios::out);

  sample_input << R"(---
initial_conditions:
  type: "ccake" # Possible values are: ["ccake", "iccing", "hdf5"]
  file: "initial_conditions/ccake_trento.dat" # Path to the initial conditions file
  t0: 0.6 # tau0 = \sqrt(t^2-z^2) in fm/c. Must match the value for initial condition file.
  dimension: 2 # (1+1)D, (2+1)D or (3+1)D simulation
  input_as_entropy: true # if false, input will be interpreted as energy density
parameters:
  dt: 0.03 # time-step of the simulation in fm/c
  h_T:  0.500000  # Smothing parameter inside kernel, in the transverse direction (in fm)
  h_eta:  0.500000  # Smothing parameter inside kernel in the longitudinal direction
  kernel_type: "quartic_spline" # Kernel type. Possible values are: ["cubic_spline"]
  energy_cutoff: 0.1 # GeV/fm^3 - Particles with less than this energy will be deleted before starting the evolution
  buffer_particles:
    enabled: true # Set to true to include buffer particles. May help stability.
    circular: true # whether to buffer with entire grid or just circular padding.
    padding_thickness: .2 # if circular buffer is true, specifies the fractional amount
                          # of padding to add beyond the point with maximum distance from origin
eos:
  type: "table" # Possible values are: ["table", "conformal"]
  path: "Houston" # If type is "table", the path to the equation of state table.
particlization:
  enabled: true # Set to true to output particlization hypersurface
  type: "fixed_T" # Possible values are: ["fixed_T"]
  T: 140.00000 # Temperature of the particlization hypersurface (in MeV)
hydro:
  baryon_charge_enabled: false # Set to false to disable baryon charge evolution
  strange_charge_enabled: false # Set to false to disable strange charge evolution
  electric_charge_enabled: false # Set to false to disable electric charge evolution
  viscous_parameters:
    shear:
      mode: "constant" # Possible values are: ["constant", "disabled"]
      constant_eta_over_s: 0.20 # Value of the constant shear viscosity
      relaxation_mode: "default" # How shear relaxation coeffiecient is computed. Possible values are: ["default"]
    bulk:
      mode: "cs2_dependent"
      constant_zeta_over_s: 0.0
      cs2_dependent_zeta_A: 1.67552
      cs2_dependent_zeta_p: 2.0
      relaxation_mode: "default"
      modulate_with_tanh: true #If true, will modulate the bulk viscosity
                               #with a tanh function below transition temperature
)";
sample_input.close();

}


int main(int argc, char* argv[]) {
 
    // Initialize Kokkos
    Kokkos::ScopeGuard guard(argc, argv);

    // Run the tests
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}