#ifndef SETTINGS_H
#define SETTINGS_H

#include <algorithm>
#include <string>
#include <vector>

using std::string;

class Settings
{
  public:

    bool using_Gubser            = false;
    bool using_Gubser_with_shear = false;
    bool using_shear             = false;

    /* DOUBLES */
    // maximum upper limit for t
    static constexpr double tend=50.02;
    
//    double step = 0.0;
    double t0 = 0.0;
    double dt = 0.0;
    double stepx = 0.0;
    double stepy = 0.0; // initial grid coordinate spacings
    double h = 0.0;
    double Freeze_Out_Temperature = 0.0;

//    double factor = 0.0;

    double efcheck = 0.0;
    double sfcheck = 0.0;



    static constexpr int VERBOSE = 5;

    string IC_type = "";        // specify initial condition type
    string IC_option = "";      // specify option for given initial condition type
    string IC_file = "";        // specify option for given initial condition type
    string EoS_type = "";       // specify equation of state type
    string EoS_option = "";     // specify specifc option for EOS
    string eta = "";            // specify the shear viscosity type to use
    string etaOption = "";      // specify necessary eta options
    string shearRelax = "";     // specify which shear relaxation to use
    string zeta = "";           // specify the bulk viscosity type to use
    string zetaOption = "";     // specify necessary zeta options
    string bulkRelax = "";      // specify which bulk relaxation to use
    string Freeze_Out_Type = "";


    vector<string> headers;


    // allows for explicitly printing tons of information about specific particles
    vector<int> particles_to_print;
    vector<bool> is_printable;

    // make sure that all chosen settings make reasonable sense
    void check_consistency()
    {
      // put any necessary consistency checks here

      return;
    }

  //  inline bool print_particle( int i ) { return is_printable[ i ]; }

};


typedef std::shared_ptr<Settings> pSettings;  // smart pointer to settings object



#endif