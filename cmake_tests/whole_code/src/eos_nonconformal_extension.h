#ifndef EOS_NONCONFORMAL_EXTENSION_H
#define EOS_NONCONFORMAL_EXTENSION_H

#include <vector>

#include "constants.h"
#include "eos_base.h"
#include "eos_extension.h"
#include "eos_header.h"

class EoS_nonconformal_extension: public EoS_base
{
  private:

    pEoS_base p_reference_EoS = nullptr;  // pointer to EoS_base object of which
                                          // this object provides an extension

  public:
    // default constructor/destructor
    EoS_nonconformal_extension(){}
    virtual ~EoS_nonconformal_extension(){}

    EoS_nonconformal_extension( pEoS_base p_reference_EoS_in,
                   const std::vector<double> & tbqs_minima_in,
                   const std::vector<double> & tbqs_maxima_in )
      { p_reference_EoS = p_reference_EoS_in;
        tbqs_minima = tbqs_minima_in; tbqs_maxima = tbqs_maxima_in;
        name = "default_with_non_conformal_extension"; }

    void get_eBSQ( double point[], double results[] )
    {
      // point: (T, muB, muQ, muS)
    }

    void get_sBSQ( double point[], double results[] )
    {
      // save original point for use below
      double point_save[4];
      for (int i = 0; i < 4; i++) point_save[i] = point[i];

      // get default EoS's grid ranges
      const auto & minima = p_reference_EoS->get_tbqs_minima_no_ext();
      const auto & maxima = p_reference_EoS->get_tbqs_maxima_no_ext();

      // project point to boundary of reference_EoS (resets point)
      eos_extension::project_to_boundary( point, minima.data(), maxima.data() );

      // set reference thermodynamics on projected point (on boundary)
      p_reference_EoS->get_sBSQ( point, results );

      // set coefficients from reference thermodynamics
      eos_extension::set_coeffs( point, results );

      // use non-conformal extension to obtain final answer (store in results)
      eos_extension::get_nonconformal_extension( point_save, point, results );
      
    }

    void get_full_thermo( double point[], double results[] )
    {

    }
    
};

#endif