#ifndef EOM_H
#define EOM_H

#include <memory>
#include <string>
#include <vector>

#include "constants.h"
#include "densities.h"
#include "hydrodynamic_info.h"
#include "settings.h"
#include "thermodynamic_info.h"

namespace ccake{
// THIS IS THE BASE CLASS
template<unsigned int D>
class EquationsOfMotion
{
  // allow derived classes to access
  protected:
    Settings * settingsPtr   = nullptr;

  public:
    EquationsOfMotion(){}
    virtual ~EquationsOfMotion(){}

    std::string name = "";                    // name associated to EoM

    void set_SettingsPtr(Settings * settingsPtr_in) { settingsPtr = settingsPtr_in; }

    // !!!require this function to be defined!!!
    virtual void evaluate_time_derivatives( hydrodynamic_info<D> & hi,
                                            thermodynamic_info & ti,
                                            densities & d_dt_specific){}
 

};

//typedef std::shared_ptr<EquationsOfMotion> pEquationsOfMotion;
}
#endif