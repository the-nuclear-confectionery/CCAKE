#ifndef EOM_BASE_H
#define EOM_BASE_H

#include <memory>
#include <string>
#include <vector>

#include "constants.h"
#include "hydrodynamic_info.h"

class EoM_base
{
  private:
    Settings * settingsPtr   = nullptr;

  public:
    EoM_base(){}
    virtual ~EoM_base(){}

    std::string name = "";                    // name associated to EoM

    void set_SettingsPtr(Settings * settingsPtr_in) { settingsPtr = settingsPtr_in; }

    // require all of these to be defined
    virtual void compute_du_dt(){}
    virtual void compute_dshv_dt(){}
    virtual void compute_detasigma_dt(){}
    virtual void compute_dBulk_dt(){}

    virtual void compute_time_derivatives(const hydrodynamic_info & h_i){}

  // allow derived classes to access and redefine
  //protected:


  //private:
    

};

typedef std::shared_ptr<EoM_base> pEoM_base;  // pointer to the base class from
                                              // which all EoMs are derived

#endif