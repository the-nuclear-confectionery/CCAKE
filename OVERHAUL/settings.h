#ifndef SETTINGS_H
#define SETTINGS_H

class Settings
{
public:

  static constexpr int VERBOSE = 5;

  //static constexpr double q=1.;
  double step;

  // maximum upper limit for t
  static constexpr double tend=50.02;

  double t0;


  double _h;
  //int _n;
  int start, end, fnum;
  int qmf;  // if==1 quantum mechanicanical corrections to the flow or added,
            // if==0 no corrections are included

  int number_part;
  double t,dt;
  double factor;
  int frzc;
  double tau, taup, taupp;
  int rk2;
  double gd2;

  int gtyp;

  int cfon;
  int cf;
  int visc; // visc=0 for ideal
            // visc=1 for bulk,
            // visc=2 for shear
            // visc=3 for bulk+shear
            // visc=4 for bulk+shear+BSQ
  double efcheck, sfcheck;
  int steps;

  int first;
  int average;
  int lowT;

  double avgetasig;
  double wfz, cs2;

  double E1, E2;

  int etaconst;
  double bvf, svf, zwidth, sTc, zTc;

  string eos_s, eos_p;
  string ebe_folder;

  int fcount, cevent;
  string eost;

  // make sure that all chosen settings make reasonable sense
  void check_consistency()
  {
    // put any necessary consistency checks here

    return;
  }

  struct Input_Parameters
{
    string IC_type; // specify initial condition type
    double h; // static SPH cutoff paramter
    double dt; // time step in fm
    double t0; // initial time in fm
    string EoS_type; // specify equation of state type
    string EoS_option; // specify specifc option for EOS
    // there should an associated EoS directory with tables
    string eta; // specificy the shear viscosity type to use
    // in transport cpefficient file
    string zeta; // specificy the bulk viscosity type to use
    // in transport cpefficient file
    double Freeze_Out_Temperature;
    string Freeze_Out_Type;
};
  Input_Parameters input_parameters;

  vector<string> headers;


};





#endif