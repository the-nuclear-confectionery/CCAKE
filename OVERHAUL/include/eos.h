#pragma once

#include <iostream>
#include <fstream>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "eos_base.h"
#include "eos_conformal_diagonal.h"
#include "eos_extension.h"
#include "eos_header.h"
#include "eos_delaunay.h"
#include "interpolatorND.h"
#include "rootfinder.h"
#include "settings.h"
#include "thermodynamic_info.h"

using std::string;

class EquationOfState
{
friend class InputOutput;

private:
    string default_eos_name = "";

public:
    ////////////////////////////////////////////////////////////////////////////
    // PUBLIC METHODS

    //Constructors:
    EquationOfState();
    EquationOfState(string quantityFile, string derivFile);

    // object to access appropriate EoS by name
    std::unordered_map<std::string, pEoS_base> chosen_EOS_map;


    ////////////////////////////////////////////////////////////////////////////
    void init();
    void init(string quantityFile, string derivFile);
//    void init_grid_ranges_only(string quantityFile, string derivFile);
    bool point_not_in_range( double setT, double setmuB, double setmuQ, double setmuS, pEoS_base peos );
    void tbqs( double setT, double setmuB, double setmuQ, double setmuS, const string & eos_name );
    void tbqs( vector<double> & tbqsIn, const string & eos_name );

    vector<double> get_thermodynamics( vector<double> & tbqsIn, const string & eos_name );


    ////////////////////////////////////////////////////////////////////////////
    //getter functions for the quantities of interest at the current tbs/tbqs
    double T()   const;   //temperature
    double muB() const;   //baryon chemical potential
    double muQ() const;   //charge chemical potential
    double muS() const;   //strangeness chemical potential

    double p()   const;   //pressure density
    double s()   const;   //entropy density
    double B()   const;   //baryon density
    double S()   const;   //strangeness density
    double Q()   const;   //charge density
    double e()   const;   //energy density
    double cs2() const;   //speed of sound
    double w()   const;   //enthalpy

    double dwds();
    double dwdB();
    double dwdS();
    double dwdQ();

    void eosin(std::string type);
    double A();

    void set_thermo( thermodynamic_info & thermo )
    {
      thermo.eos_name = get_current_eos_name();

      thermo.T    = T();
      thermo.muB  = muB();
      thermo.muS  = muS();
      thermo.muQ  = muQ();

      thermo.p    = p();
      thermo.s    = s();
      thermo.rhoB = B();
      thermo.rhoS = S();
      thermo.rhoQ = Q();
      thermo.e    = e();
      thermo.w    = w();
      thermo.A    = A();
      thermo.cs2  = cs2();
      thermo.dwds = dwds();
      thermo.dwdB = dwdB();
      thermo.dwdS = dwdS();
      thermo.dwdQ = dwdQ();
    }

    // call these functions using default EoS if none is specified
    double efreeze(double TFO){return efreeze(TFO, default_eos_name);}
    double sfreeze(double TFO){return sfreeze(TFO, default_eos_name);}
    double cs2out(double Tt, double muBin, double muQin, double muSin)
           {return cs2out(Tt, muBin, muQin, muSin, default_eos_name);}
    double cs2out(double Tt){return cs2out(Tt, default_eos_name);}
    double wfz(double Tt, double muBin, double muQin, double muSin)
           {return wfz(Tt, muBin, muQin, muSin, default_eos_name);}
    double wfz(double Tt){return wfz(Tt, default_eos_name);}
    double s_terms_T(double Tt){return s_terms_T(Tt, default_eos_name);}

    double efreeze(double TFO, const string & eos_name);
    double sfreeze(double TFO, const string & eos_name);
    double cs2out(double Tt, double muBin, double muQin, double muSin,
                  const string & eos_name);
    double cs2out(double Tt, const string & eos_name);
    double wfz(double Tt, double muBin, double muQin, double muSin,
                  const string & eos_name);
    double wfz(double Tt, const string & eos_name);
    double s_terms_T(double Tt, const string & eos_name);

    double s_terms_T( double Tt, double muBin, double muQin, double muSin,
                      const string & eos_name );

    void evaluate_thermodynamics( pEoS_base peos );

    bool update_s(double sin, double Bin, double Sin, double Qin);
    bool update_s(double sin);
    double s_out(double ein, double Bin, double Sin, double Qin, bool & solution_found);
    double s_out(double ein, bool & solution_found);

    void set_SettingsPtr( Settings * settingsPtr_in );

    vector<pEoS_base> chosen_EOSs;     // the vector of EoSs to use, in order

    // use this specifically to make educated guesses for seed values, etc.
//    EoS_conformal_diagonal conformal_diagonal_EoS;

private:

    ////////////////////////////////////////////////////////////////////////////
    // PRIVATE MEMBERS

    Settings * settingsPtr = nullptr;

    //the current position in (T, muB, muQ, muS) initialized by tbqs()
    vector<double> tbqsPosition;

    // string to hold input filenames and EoS table
    string quantity_file = "";
    string deriv_file    = "";
    string current_eos_name = "";

    double pVal          = 0.0;
    double entrVal       = 0.0;
    double BVal          = 0.0;
    double SVal          = 0.0;
    double QVal          = 0.0;
    double eVal          = 0.0;
    double cs2Val        = 0.0;

    double db2           = 0.0;
    double ds2           = 0.0;
    double dq2           = 0.0;
    double dt2           = 0.0;
    double dbdq          = 0.0;
    double dbds          = 0.0;
    double dsdq          = 0.0;
    double dtdb          = 0.0;
    double dtds          = 0.0;
    double dtdq          = 0.0; //second derivative of pressure wrt i and j 
                                //where didj =: (d^2p)/(didj) or di2 = (d^2p)/((di)^2)


    ////////////////////////////////////////////////////////////////////////////
    // PRIVATE ROUTINES FOR SETTING POINT IN PHASE DIAGRAM FOR GIVEN EOS
    void tbqs(double setT, double setmuB, double setmuQ, double setmuS, pEoS_base peos);
    void tbqs( vector<double> & tbqsIn, pEoS_base peos )
          { tbqs( tbqsIn[0], tbqsIn[1], tbqsIn[2], tbqsIn[3], peos ); }



    ////////////////////////////////////////////////////////////////////////////
    // ROUTINES NEEDED FOR COMPUTING THERMODYNAMIC DERIVATIVES
    double dentr_dt();
    double dentr_dmub();
    double dentr_dmuq();
    double dentr_dmus();
    double db_dt();
    double db_dmub();
    double db_dmuq();
    double db_dmus();
    double ds_dt();
    double ds_dmub();
    double ds_dmuq();
    double ds_dmus();
    double dq_dt();
    double dq_dmub();
    double dq_dmuq();
    double dq_dmus();
    double calc_term_1();
    double calc_term_2(string i_char);
    double calc_term_3(string i_char);
    double calc_term_4(string j_char, string i_char);
    double deriv_mult_aTm_1b(gsl_vector* a, gsl_matrix* m, gsl_vector* b);

    
    ////////////////////////////////////////////////////////////////////////////
    // ROUTINES FOR DEBUGGING (uncomment to call)
    //void check_EoS_derivatives();
    //void get_toy_thermo(double point[], double thermodynamics[]);


    ////////////////////////////////////////////////////////////////////////////
    // MISCELLANEOUS PRIVATE ROUTINES
    bool find_root_with_seed( const string & e_or_s_mode, double e_or_s_in,
                              double Bin, double Sin, double Qin,
                              pEoS_base this_eos, vector<double> & result );

    bool delaunay_update_s(double sin, double Bin, double Sin, double Qin);
    bool rootfinder_update_s(double sin, double Bin, double Sin, double Qin);
    double delaunay_s_out(double ein, double Bin, double Sin, double Qin, bool & solution_found);
    double rootfinder_s_out(double ein, double Bin, double Sin, double Qin, bool & solution_found);

    // need these to be static to initialize std::function<...> objects
    static void get_eBSQ_densities_from_interpolator( double point[], double densities[] );
    static void get_sBSQ_densities_from_interpolator( double point[], double densities[] );


    ////////////////////////////////////////////////////////////////////////////
    // MEMBERS AND ROUTINES TO FIND (T,muX) COORDINATES OF (e,rhoX) POINT
    // - for using the root-finding functionality
    Rootfinder rootfinder;
    // - for using a Delaunay interpolation
    eos_delaunay e_delaunay;
    eos_delaunay entr_delaunay;

public:
  //bool using_conformal_as_fallback() { return use_conformal_as_fallback; }
  string get_current_eos_name() { return current_eos_name; }
  void run_closure_test();

  //============================================================================
  // function to check whether conformal-diagonal EoS instance
  // yields a real solution for this (e,B,S,Q) combination
  bool eBSQ_has_solution_in_conformal_diagonal(
        double e_in, double B_in, double S_in, double Q_in )
  {
    // select the currently defined conformal-diagonal EoS instance
    auto conformal_diagonal_EoS
            = std::dynamic_pointer_cast<EoS_conformal_diagonal>
                        ( chosen_EOS_map["conformal_diagonal"] );
    return conformal_diagonal_EoS->eBSQ_has_solution( e_in, B_in, S_in, Q_in );
  }

};
