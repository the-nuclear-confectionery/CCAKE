#ifndef EOM_H
#define EOM_H

#include <memory>
#include <string>
#include <vector>

#include "constants.h"
#include "hydrodynamic_info.h"
#include "settings.h"

// THIS IS THE BASE CLASS
class EquationsOfMotion
{
  // allow derived classes to access and redefine
  protected:
    Settings * settingsPtr   = nullptr;

  public:
    EquationsOfMotion(){}
    virtual ~EquationsOfMotion(){}

    std::string name = "";                    // name associated to EoM

    void set_SettingsPtr(Settings * settingsPtr_in) { settingsPtr = settingsPtr_in; }

    // require all of these to be defined
    virtual void compute_du_dt(){}
    virtual void compute_dshv_dt(){}
    virtual void compute_detasigma_dt(){}
    virtual void compute_dBulk_dt(){}

    virtual void evaluate_time_derivatives(const hydrodynamic_info & h_i)
{cout << "t=: In " << __LINE__ << endl;}

  //private:
    

};

typedef std::shared_ptr<EquationsOfMotion> pEquationsOfMotion;

#endif