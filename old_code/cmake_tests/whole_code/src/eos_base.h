#ifndef EOS_BASE_H
#define EOS_BASE_H

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "constants.h"
#include "eos_header.h"

class EoS_base
{
  public:
    EoS_base(){}
    virtual ~EoS_base(){}

    // functionals every EoS needs
    std::function<void(double[], double[])> eBSQ
              = [this](double a[], double b[]) { this->get_eBSQ(a, b); };
    std::function<void(double[], double[])> sBSQ
              = [this](double a[], double b[]) { this->get_sBSQ(a, b); };
    std::function<void(double[], double[])> full_thermo
              = [this](double a[], double b[]) { this->get_full_thermo(a, b); };

    //std::function<void(int,int)> f = [this](int a, int b)
    //        { this->doSomethingArgs(a, b); };

    std::vector<double> tbqs_minima, tbqs_maxima; // ranges over which EoS is defined
                                                  // (including possible extension)

    std::string name = "";                    // name associated to EoS

    virtual void get_eBSQ( double point[], double results[] ){}
    virtual void get_sBSQ( double point[], double results[] ){}
    virtual void get_full_thermo( double point[], double results[] ){}

    const std::vector<double> & get_tbqs_minima_no_ext()
                                const { return tbqs_minima_no_ext; }
    const std::vector<double> & get_tbqs_maxima_no_ext()
                                const { return tbqs_maxima_no_ext; }

  // allow derived classes to access and redefine
  protected:
    std::vector<double> tbqs_minima_no_ext, tbqs_maxima_no_ext;
                        // ranges over which EoS is defined
                        // (excluding possible extension,
                        // may be the same as above)

};

typedef std::shared_ptr<EoS_base> pEoS_base;  // pointer to the base class from
                                              // which all EoSs are derived

#endif