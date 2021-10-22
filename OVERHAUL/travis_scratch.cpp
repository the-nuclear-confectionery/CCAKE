///////////INPUT_OUTPUT.CPP BEGIN/////////////

#ifndef INPUT_OUTPUT_CPP

void load_settings_file( string path_to_settings_file )
{
    string Param_file = path_to_settings_file+"Input_Parameters.inp"
    ifstream infile( Param_file.c_str() );
    if (infile.is_open())
    {
        string line;
        string ignore
        vector <string> param;
        int param_line = 0
        while ( getline (infile, line) )
        {
            istringstream iss(line);
            iss >> ignore >> param[param_line];
        }

        input_parameter.IC_tpye = param[0]
        input_parameter.h = stod(param[1])
        input_parameter.dt = stod(param[2])
        input_parameter.t0 = stod(param[3])
        input_parameter.EoS = param[4]
        input_parameter.eta = param[5]
        input_parameter.zeta = param[6]
        input_parameter.Freeze_Out_Temperature = stod(param[7])
        input_parameter.Freeze_Out_Type = param[8]

        infile.close();
    }
}

///////////INPUT_OUTPUT.CPP END//////////////


///////////BSQHYDRO.CPP BEGIN/////////////

#include "BSQHydro.h"
#include "input_output.h"


// Constructors and destructors.
  BSQHydro::BSQHydro(){}
 ~BSQHydro::BSQHydro(){}


void BSQHydro::load_settings_file( string path_to_settings_file ){}



void BSQHydro::set_results_directory( string path_to_results_directory ){}



void BSQHydro::read_in_initial_conditions(){}



void BSQHydro::initialize_hydrodynamics(){}


void BSQHydro::run(){}


void BSQHydro::find_freeze_out_surface(){}


void BSQHydro::print_results(){}
///////////BSQHYDRO.CPP END//////////////





///////////BSQHYDRO.H BEGIN/////////////

#ifndef BSQHYDRO_H
#define BSQHYDRO_H

using std::string;

class BSQHydro
{

public:

  BSQHydro();
  ~BSQHydro();

  void load_settings_file( string path_to_settings_file );
  void set_results_directory( string path_to_results_directory );
  void read_in_initial_conditions();
  void initialize_hydrodynamics();
  void run();
  void find_freeze_out_surface();
  void print_results();

  struct Input_Parameters
  {
    string IC_type; // specify initial condition type
    double h; // static SPH cutoff paramter
    double dt; // time step in fm
    double t0; // initial time in fm
    string EoS; // specify equation of state type
    // there should an associated EoS directory with tables
    string eta; // specificy the shear viscosity type to use
    // in transport cpefficient file
    string zeta; // specificy the bulk viscosity type to use
    // in transport cpefficient file
  }

  struct Initial_Conditions
  {
      vector<vector string> headers;
      vector<vector double> density_grid;
  }



private:

  string input_directory
  string output_directory;

  // equation of state information
  EquationOfState eos;

  // input/output
  InputOutput io;

  // hold freeze-out surface
  FreezeOutSurface freeze_out_surface;

  // the current state of the hydrodynamic simulation
  SystemState system;

}

#endif

///////////BSQHYDRO.H END//////////////
