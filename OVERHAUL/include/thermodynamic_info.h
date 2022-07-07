#ifndef THERMODYNAMIC_INFO_H
#define THERMODYNAMIC_INFO_H


// thermodynamic quantities (with default initialization)
struct thermodynamic_info
{
  // (T,mu_i) coordinates depend on which EoS was used!
  string eos_name = "";

  double T    = 0.0, muB  = 0.0, muS  = 0.0, muQ  = 0.0;
  double e    = 0.0, s    = 0.0, rhoB = 0.0, rhoS = 0.0, rhoQ = 0.0,
         p    = 0.0, cs2  = 0.0, w    = 0.0, A    = 0.0;
  double dwds = 0.0, dwdB = 0.0, dwdS = 0.0, dwdQ = 0.0;

};


#endif