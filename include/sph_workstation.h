#ifndef SPH_WORKSTATION_H
#define SPH_WORKSTATION_H

#include "eom.h"
#include "eom_default.h"
#include "evolver.h"
#include "formatted_output.h"
#include "freeze_out.h"
#include "kernel.h"
#include "settings.h"
#include "stopwatch.h"
#include "system_state.h"
#include "transport_coefficients.h"

class SPHWorkstation
{
friend class BSQHydro;
friend class Evolver;
friend class InputOutput;

private:
  
  static constexpr int    VERBOSE        = 3;
  static constexpr double TOLERANCE      = 0.0;
  static constexpr bool   REGULATE_LOW_T = false;

  Settings        * settingsPtr          = nullptr;
  SystemState     * systemPtr            = nullptr;

  // equations of motion
  pEquationsOfMotion pEoM;

  // evolver
  Evolver evolver;

  // equation of state
  EquationOfState eos;

  // transport coefficients
  TransportCoefficients tc;

  // freeze out
  FreezeOut fo;


public:

  //============================================================================
  // default constructor/destructor
  SPHWorkstation(){}
  ~SPHWorkstation(){}

  //============================================================================
  // initialize pointers
  void set_SystemStatePtr( SystemState * systemPtr_in );
  void set_SettingsPtr( Settings * settingsPtr_in );

  //============================================================================
  // initialize workstation (now includes eos initialization)
  void initialize()
  {
    //----------------------------------------
    // set up equation of motion object
    pEoM = std::make_shared<EoM_default>();
    pEoM->set_SettingsPtr( settingsPtr );

    //----------------------------------------
    // set up equation of state
    eos.set_SettingsPtr( settingsPtr );
    eos.init();

    //----------------------------------------
    // set up transport coefficients
    tc.set_SettingsPtr( settingsPtr );
    //tc.initialize();  // this just uses the default transport coefficient modes
    if ( settingsPtr->using_Gubser )
      tc.initialize( "Gubser" );
    else
      tc.initialize( settingsPtr->etaMode, settingsPtr->shearRelaxMode,
                     settingsPtr->zetaMode, settingsPtr->bulkRelaxMode );

    //----------------------------------------
    // set up freeze out (constant energy density
    fo.set_SettingsPtr( settingsPtr );
    fo.set_SystemStatePtr( systemPtr );
    systemPtr->efcheck = eos.efreeze(settingsPtr->Freeze_Out_Temperature);
    systemPtr->sfcheck = eos.sfreeze(settingsPtr->Freeze_Out_Temperature);
    fo.initialize( systemPtr->efcheck );

    //----------------------------------------
    // set up RK evolver
    evolver.set_SettingsPtr( settingsPtr );
    evolver.set_SystemStatePtr( systemPtr );
  }

  //============================================================================
  // routines for resetting quantities
  void reset_linklist() { systemPtr->linklist.reset(); }
  void reset_pi_tensor();

  void process_initial_conditions();
  void initialize_entropy_and_charge_densities();
  void initial_smoothing();
  void add_buffer(double default_e);

  //============================================================================
  // smoothing
  void smooth_fields( Particle & pa );
  void smooth_gradients( Particle & pa, double tin );

  void get_time_derivatives();


  //============================================================================
  // functions to apply action to all particles
  //-------------------------------------------
  void smooth_all_particle_fields()
        { Stopwatch sw;
          sw.Start();
          for ( auto & p : systemPtr->particles )
            smooth_fields(p);
          sw.Stop();
          formatted_output::update("finished smoothing particle fields in "
                                    + to_string(sw.printTime()) + " s.");
        }

  //-------------------------------------------
  void smooth_all_particle_gradients()
        {
          Stopwatch sw;
          sw.Start();
          for ( auto & p : systemPtr->particles )
            smooth_gradients( p, systemPtr->t );
          sw.Stop();
          formatted_output::update("finished smoothing particle gradients in "
                                    + to_string(sw.printTime()) + " s.");
        }

  //-------------------------------------------
  void update_all_particle_thermodynamics()
        {
          Stopwatch sw;
          sw.Start();
          for ( auto & p : systemPtr->particles )
            set_phase_diagram_point( p );
          sw.Stop();
          formatted_output::update("got particle thermodynamics in "
                                    + to_string(sw.printTime()) + " s."); }

  //-------------------------------------------
  void update_all_particle_viscosities()
        { for ( auto & p : systemPtr->particles )
            set_transport_coefficients( p ); }

  //-------------------------------------------
  void evaluate_all_particle_time_derivatives()
        { 
          Stopwatch sw;
          sw.Start();
          for ( auto & p : systemPtr->particles )
          {
            p.hydro.ID = p.ID;
            p.hydro.t  = systemPtr->t;
            pEoM->evaluate_time_derivatives( p.hydro, p.thermo, p.d_dt_spec );
          }
          sw.Stop();
          formatted_output::update("set particle time derivatives in "
                                    + to_string(sw.printTime()) + " s."); }


  //-------------------------------------------
  void freeze_out_particles();


  //============================================================================
  // routines to edit particles directly
  double locate_phase_diagram_point_eBSQ(
          Particle & p, double e_In,
          double rhoB_In, double rhoS_In, double rhoQ_In );
  double locate_phase_diagram_point_eBSQ(
          Particle & p, double e_In );
  void locate_phase_diagram_point_sBSQ(
          Particle & p, double s_In,
          double rhoB_In, double rhoS_In, double rhoQ_In );
  void locate_phase_diagram_point_sBSQ(
          Particle & p, double s_In );

  void set_phase_diagram_point( Particle & p );
  void set_transport_coefficients( Particle & p );

  void set_bulk_Pi();
  void set_dcs2_dt();//intialize dcs2_dt at first time step
  // misc. routine
  double gradPressure_weight(const int a, const int b);
  double gradEnergy_weight(const int a, const int b);

  //============================================================================
  // what it says on the label
  void advance_timestep( double dt, int rk_order )
  {
    Stopwatch sw;
    sw.Start();

    // turn on freeze-out flag initially
    systemPtr->do_freeze_out = true;

    // use evolver to actually do RK evolution
    // (pass workstation's own time derivatives function as lambda)
    evolver.advance_timestep( dt, rk_order,
                              [this]{ this->get_time_derivatives(); } );

    // set number of particles which have frozen out
    systemPtr->number_part = systemPtr->get_frozen_out_count();

    // keep track of how many timesteps have elapsed
    systemPtr->number_of_elapsed_timesteps++;

    sw.Stop();
    formatted_output::report("finished timestep in "
                              + to_string(sw.printTime()) + " s");

    if ( settingsPtr->max_number_of_timesteps >= 0
          && systemPtr->number_of_elapsed_timesteps
              > settingsPtr->max_number_of_timesteps )
      exit(-1);

    return;
  }

  //============================================================================
  // decide whether to continue evolving
  bool continue_evolution()
  {
    return ( systemPtr->t < settingsPtr->tend )
            && ( systemPtr->number_part < systemPtr->n() );
  }

};



#endif