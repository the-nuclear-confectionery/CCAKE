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

        input_parameters.IC_tpye = param[0];
        input_parameters.h = stod(param[1]);
        input_parameters.dt = stod(param[2]);
        input_parameters.t0 = stod(param[3]);
        input_parameters.EoS_type = param[4];
        input_parameters.EoS_option = param[5];
        input_parameters.eta = param[6]
        input_parameters.zeta = param[7]
        input_parameters.Freeze_Out_Temperature = stod(param[8])
        input_parameters.Freeze_Out_Type = param[9]

        infile.close();
    }
}

void set_EoS_type()
{
  EoS_type = input_parameters.EoS_type;
  string EoS_files_location = 'EoS/' + EoS_type;
  string densities = EoS_files_location + '/densities.dat';
  string derivatives = EoS_files_location + '/derivatives.dat';
  string EoS_option = input_paramters.EoS_option;

  switch(EoS_spec)
  {
    Case default   :
      cout << "Running default EoS option for " << EoS_type << endl;
  }

  eos.innit(desnities,derivatives);
  return
}



void InputOutput::read_in_initial_conditions()
{
  string initial_condition_type = input_parameters.IC_type;
  int total_header_lines;
  string IC_file = 'All_Initial_Conds/'
  switch(initial_condition_type)
  {
    case 'ICCING' :
      cout << "Reading in ICCING initial conditions!" << endl;
      IC_file = IC_file+'Iccing_conditions.dat' // need to change ic0.dat
      total_header_lines = 1;
    case default :
      cout << "Selected initial condition type not supported."
  }

  ifstream infile(IC_file.c_str());
  cout << "Initial conditions file: " << IC_file << endl;
  if (infile.is_open())
  {
    string line;
    int count_header_lines = 0;
    int count_file_lines = 0;
    double x,y,e,rhoB,rhoS,rhoQ;
  while (getline (infile, line))
          {
            if(count_header_lines < total_header_lines)
            {
              initial_conditions.headers.pushback(line)
              count_header_lines++;
            }
            else
            {
              istringstream iss(line);
              iss >> x >> y >> e >> rhoB >> rhoS >> rhoQ;
              vector<double> fields({x,y,e,rhoB,rhoS,rhoQ})
              initial_conditions.density_grid.pushback(fields)
          }
        }
  }
  else
  {

    cout << "Can't open " << IC_file << endl;
    exit(1);
  }
  infile.close();
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



void BSQHydro::trim_initical_conditions()
{
  vector<vector<double>> trimmed_grid;
  int cells_before_trim = initial_conditions.density_grid.size();
  for(int i=0; i<cells_before_trim; i++)
  {
    double e = initial_conditions.density_grid[i][2];
    double rhoB = initial_conditions.density_grid[i][3];
    double rhoS = initial_conditions.density_grid[i][4];
    double rhoQ = initial_conditions.density_grid[i][5];

    eos.sout(e,rhoB,rhoS,rhoQ);
    if (eos.T() > input_parameters.Freeze_Out_Temperature)
    {
      trimmed_grid.pushback(initial_conditions.density_grid[i]);
    }
  }
  initial_conditions.density_grid = trimmed_grid;
}



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


  void trim_initical_conditions()

}

#endif

///////////BSQHYDRO.H END//////////////
