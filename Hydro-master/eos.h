#pragma once

#include "../../splinter/include/bspline.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>

#include <fstream>
#include <string>

using std::string;

using namespace SPLINTER;

class eos {

public:

    //Constructors:
    //always either use the eos(string, string, int) constructor or the default constructor followed by the init function
    //both of those actions will accomplish the same thing

    //quantityFile must be in dimensionless quantities and must be formatted "T  muB  muQ  muS  p  s  B  S  Q  e  cs2"
    //derivFile should be formatted "T  muB  muQ  muS  d2p/dB2  d2p/dQ2  d2p/dS2  d2p/dBdQ  d2p/dBdS d2p/dQdS  d2p/dTdB  d2p/dTdQ  d2p/dTdS  d2p/dT2"
    eos(string quantityFile, string derivFile, int degree = 4);

    eos();
    void init(string quantityFile, string derivFile, int degree = 4);


    //initializes the position in the grid to (setT,setmuB,setmuQ,setmuS)
    //Once called, the splines will stay initialized at this point until the function is called again
    void tbqs(double setT, double setmuB, double setmuQ, double setmuS);

    //getter functions for the quantities of interest at the current tbs/tbqs
    double p();     //pressure density
    double s();     //entropy density
    double B();     //baryon density
    double S();     //strangeness density
    double Q();     //charge density
    double e();     //energy density
    double cs2();   //speed of sound
    double w();     //enthalpy
    double T();     //temperature
    double muB();   //baryon chemical potential
    double muQ();   //charge chemical potential
    double muS();   //strangeness chemical potential

    double dwds();
    double dwdB();  //enthalpy derivatives **These still have not been checked**
    double dwdS();
    double dwdQ();

    double cs2out(double Tt, double muBin, double muQin, double muSin); //return cs2 given T and mu's - identical to calling cs2() after initializing position using tbqs()
    double cs2out(double Tt); //return cs2 given T and mu's = 0 - identical to calling cs2() after initializing position using tbqs()
    double wfz(double Tt, double muBin, double muQin, double muSin);    // return enthalpy for T and mu's - identical to calling w() after initializing position using tbqs()
    double wfz(double Tt);    // return enthalpy for T and mu's = 0 - identical to calling w() after initializing position using tbqs()
    double s_terms_T(double Tt);                              //returns entropy at a given T for muB = muS = muQ = 0


    //Rootfinder functions:

    //finds and initializes position in T and mu's with given s,BSQ - rootfinder
    //Updates the t and mu position and returns 1 if found, returns 0 if failed
    //Runs using the initial guess of the previous tbqs call (or previous result of rootfinder run)
    //If the initial guess does not succeed, rootfinder perturbs the initial guess in every T, mu direction by 20% of the initial guess
    //If any of the perturbed solutions succeed, the update_s will return 1, otherwise returns 0
    //If rootfinder fails, the tbqs position defaults to the previous value
    bool update_s(double sin, double Bin, double Sin, double Qin);
    // for backwards compatibility
    bool update_s(double sin);

    //finds and initializes position in T and mu's with given e,BSQ - rootfinder
    //Updates the t and mu position and returns entropy if found, returns -1 if failed
    //Runs using the initial guess of the previous tbqs call (or previous result of rootfinder run)
    //If the initial guess does not succeed, rootfinder perturbs the initial guess in every T, mu direction by 20% of the initial guess
    //If any of the perturbed solutions succeed, the s_out will return the entropy, otherwise returns -1
    //If rootfinder fails, the tbqs position defaults to the previous value
    double s_out(double ein, double Bin, double Sin, double Qin);
    // for backwards compatibility
    double s_out(double ein);



    // //**These functions are placeholders for the sake of compilation. They do nothing and should not be used**
    void eosin(std::string type);
    double A();

    double efreeze();
    double sfreeze();

private:

    //see READMEeos file for more info on the following parameters
    //tolerance for rootfinder function used in update_s() and s_out() **THIS CAN BE CHANGED DEPENDING ON EOS
    double TOLERANCE = 1e-11;
    //max number of steps taken by rootfinder before failure is declared **THIS CAN BE CHANGED DEPENDING ON EOS
    size_t STEPS = 1000;
    //Rootfinding method used **THIS CAN BE CHANGED DEPENDING ON EOS
    const gsl_multiroot_fsolver_type *TYPE = gsl_multiroot_fsolver_hybrids;

    //interpolations of each quantity
    BSpline pSpline, entrSpline, bSpline, sSpline, qSpline, eSpline, cs2Spline;
    BSpline db2Spline, dq2Spline, ds2Spline, dt2Spline, dbdqSpline, dbdsSpline, dtdbSpline, dqdsSpline, dtdqSpline, dtdsSpline;

    //value of each quantity at the current tbsPosition
    double pVal, entrVal, BVal, SVal, QVal, eVal, cs2Val;
    double db2, ds2, dq2, dt2, dbdq, dbds, dsdq, dtdb, dtds, dtdq; //second derivative of pressure wrt i and j where didj =: (d^2p)/(didj) or di2 = (d^2p)/((di)^2)

    double maxMuB, minMuB, maxMuQ, minMuQ, maxMuS, minMuS, maxT, minT; //EOS range used for rootfinder checks

    double dentr_dt();
    double dentr_dmub();
    double dentr_dmuq();
    double dentr_dmus();
    double db_dt();
    double db_dmub();
    double db_dmuq();
    double db_dmus();
    double ds_dt();          //intermediate derivative calculation functions
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
    double Aout;
    double Aideal();
    double Atable();
    //the current position in (T, muB, muQ, muS) initialized by tbqs()
    DenseVector tbqsPosition;

    //finds t,mu for a given e and BSQ or s and BSQ
    //returns 1 if the point was found. Returns 0 if failed
    bool rootfinder4D(double e_or_s_Given, int e_or_s_mode, double rhoBGiven, double rhoSGiven, double rhoQgiven, double error, size_t steps);
    //finds lines of constant e, p, entropy, or gibbs in the eos space
    //returns 1 if the point was found. Returns 0 if failed
    //quant type can be "e","p","entr","gibbs"
    //whichIndepVar chooses the variable dimension along the line. 1 for varying T, 2 for varyng muB, 3 for varying muQ, 4 for varying muS
    bool quant_rootfinder4D(double quantGiven, string quantType, int whichIndepVar, double error, size_t steps);
};
