#ifndef TRANSPORT_COEFFICIENTS_H
#define TRANSPORT_COEFFICIENTS_H
#include <string>
#include <vector>
#include <functional>
#include <math.h>
#include "settings.h"
#include "kernel.h"
#include "particle.h"
#include "new_matrix.h"
#include "eos.h"

using std::string;
using std::vector;

class TransportCoeficients
{

  public:
    TransportCoeficients(){};
    ~TransportCoeficients(){};

    double getEta();
    double getZeta();
    // Matrix getKappa();
    double getTauShear();
    double getTauBulk();
    // Matrix getTauDiffusive();
    void set_EquationOfStatePtr( EquationOfState * eosPtr_in );


  private:
    string etaTable, zetaTable; //read in path and load hear or should I/O
    //load in directly?? Similar qeustion for EOS..

    // struct constEta {
    //   // Constant eta functor, constructor takes in constant
    //   // etaT/w as its paramter, functor returns eta
    //   constEta(double inEta) : eta_T_OV_w(inEta) {}
    //   double operator()() const { return
    //    eta_T_OV_w_in*(eosPtr->w()/eosPtr->T()); }
    //   private:
    //     double eta_T_OV_w_in;
    // };
    // struct funcEta {
    //   // etaT/w as a function, written in the current way this
    //   // function is always formally of (T,muB,muS,muQ), but
    //   // there's no reason the function you insert must have
    //   // these depencies. Insert any function as long as it
    //   // returns a double and takes in 4 doubles
    //   // (maybe convert input to array? set the array in a
    //   // specefic order, write functions to pull as needed from array..)
    //   funcEta(function<double(double,double,double,double) func):
    //     eta_T_OV_w(func) {}
    //   double operator()() const { return 
    //   eta_T_OV_w(eosPtr->T(),eosPtr->muB(),eosPtr->muS(),eosPtr->muQ());}
    //   private:
    //     function<double(double,double,double,double) eta_T_OV_w; 
    // };
    // // Could maybe collapse above structs into one? but then the constant eta
    // // case would evaluate T,mu_i for every particle for no reason..
    // // Could also make separate structs for mu_i=0 vs mu_i =/= 0..
    // // more structs = more clutter = less evaluations
    // // less structs = less clutter = more evaluations
    // // either way it's modular and assigns the correct function at runtime
    // // without having to check every time

    // // in the case of using arrays, all of these functions would be
    // // written to take in vector<double>, then the functions would
    // // only access the elemenets they need
    // double JakiParam(double t, double mub, double mus, double muq);
    // double LinearMusParam(double t, double mub, double mus, double muq);
    // double InterpolantWrapper(double t, double mub, double mus, double muq);

    double constEta();
    double JakiParam();
    double LinearMusParam();
    double InterpolantWrapper();
    function<double()> eta;
    
    double tauShearGubser();
    double tauShearMinval();
    function<double()> tauShear;

    double zetaAdSParam();
    function<double()> zeta;

    double tauBulkDNMR_LeadingMass();
    function<double()> tauBulk;


    EquationOfState * eosPtr = nullptr;

    









}



#endif
