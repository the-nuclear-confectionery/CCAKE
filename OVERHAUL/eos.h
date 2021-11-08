#pragma once

#include <iostream>
#include <vector>
#include <fstream>
#include <functional>
#include <string>

//#include "read_in_hdf/read_in_hdf.h"
//#include "Stopwatch.h"
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "eos_delaunay/eos_delaunay.h"
#include "interpolatorND/interpolatorND.h"
#include "rootfinder.h"

using std::string;

class EquationOfState
{
friend class InputOutput;

public:
    ////////////////////////////////////////////////////////////////////////////
    // PUBLIC METHODS

    //Constructors:
    EquationOfState();
    EquationOfState(string quantityFile, string derivFile);

    void init();
    void init(string quantityFile, string derivFile);
    void init_grid_ranges_only(string quantityFile, string derivFile);
    void tbqs(double setT, double setmuB, double setmuQ, double setmuS);
    void tbqs( vector<double> & tbqsIn );

    //getter functions for the quantities of interest at the current tbs/tbqs
    double T()   const;     //temperature
    double muB() const;   //baryon chemical potential
    double muQ() const;   //charge chemical potential
    double muS() const;   //strangeness chemical potential

    double p()   const;     //pressure density
    double s()   const;     //entropy density
    double B()   const;     //baryon density
    double S()   const;     //strangeness density
    double Q()   const;     //charge density
    double e()   const;     //energy density
    double cs2() const;   //speed of sound
    double w()   const;     //enthalpy

    double dwds();
    double dwdB();  //enthalpy derivatives **These still have not been checked**
    double dwdS();
    double dwdQ();

    void eosin(std::string type);
    double A();

    double efreeze(double TFO);
    double sfreeze(double TFO);

    double cs2out(double Tt, double muBin, double muQin, double muSin);
    double cs2out(double Tt);
    double wfz(double Tt, double muBin, double muQin, double muSin);
    double wfz(double Tt);
    double s_terms_T(double Tt); 

    void evaluate_thermodynamics();

    bool update_s(double sin, double Bin, double Sin, double Qin);
    bool update_s(double sin);
    double s_out(double ein, double Bin, double Sin, double Qin);
    double s_out(double ein);

    void set_eBSQ_functional( std::function<void(double[], double[])> fIn ) { eBSQ_functional = fIn; }
    void set_sBSQ_functional( std::function<void(double[], double[])> fIn ) { sBSQ_functional = fIn; }

private:

    ////////////////////////////////////////////////////////////////////////////
    // PRIVATE MEMBERS

    //bool use_delaunay = false, use_rootfinder = true;
    //const int VERBOSE = 1;
    static constexpr bool use_rootfinder                      = true;
    static constexpr bool use_delaunay                        = !use_rootfinder;
    static constexpr bool use_static_C_library                = true;
    static constexpr bool accept_nearest_neighbor             = false;
    static constexpr bool discard_unsolvable_charge_densities = false;

    static constexpr size_t STEPS     = 1000000;
    static constexpr int VERBOSE      = 0;
    static constexpr double TOLERANCE = 1e-12;

    //the current position in (T, muB, muQ, muS) initialized by tbqs()
    vector<double> tbqsPosition;

    // string to hold input filenames and EoS table
    string quantity_file = "";
    string deriv_file    = "";
    string equation_of_state_table_filename = "";
    InterpolatorND<4> equation_of_state_table;


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

    double maxMuB        = 0.0;
    double minMuB        = 0.0;
    double maxMuQ        = 0.0;
    double minMuQ        = 0.0;
    double maxMuS        = 0.0;
    double minMuS        = 0.0;
    double maxT          = 0.0;
    double minT          = 0.0; //EOS range used for rootfinder checks

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
    bool delaunay_update_s(double sin, double Bin, double Sin, double Qin);
    bool rootfinder_update_s(double sin, double Bin, double Sin, double Qin);
    double delaunay_s_out(double ein, double Bin, double Sin, double Qin);
    double rootfinder_s_out(double ein, double Bin, double Sin, double Qin);

    // need these to be static to initialize std::function<...> objects
    void get_eBSQ_densities_from_interpolator( double point[], double densities[] );
    void get_sBSQ_densities_from_interpolator( double point[], double densities[] );


    ////////////////////////////////////////////////////////////////////////////
    // MEMBERS AND ROUTINES TO FIND (T,muX) COORDINATES OF (e,rhoX) POINT
    // - for using the root-finding functionality
    Rootfinder rootfinder;
    std::function<void(double[], double[])> eBSQ_functional;
    std::function<void(double[], double[])> sBSQ_functional;
    // - for using a Delaunay interpolation
    eos_delaunay e_delaunay;
    eos_delaunay entr_delaunay;

};
