#ifndef EOM_DEFAULT_H
#define EOM_DEFAULT_H

#include "particle.h"
#include "utilities.h"
#include "densities.h"
#include "eom.h"
#include "hydrodynamic_info.h"
#include "matrix.h"
#include "settings.h"
#include "thermodynamic_info.h"
#include "vector.h"

#include <Cabana_Core.hpp>

namespace ccake{
template<unsigned int D>
class EoM_default
{
  private:
    static constexpr int VERBOSE = 3;
    std::shared_ptr<Settings> settingsPtr;
    Matrix<double,D,D> Imat;

  public:
    EoM_default() = delete;
    EoM_default( std::shared_ptr<Settings> settingsPtr_in ): settingsPtr(settingsPtr_in) { Imat.identity(); };
    virtual ~EoM_default(){};

    KOKKOS_FUNCTION
    static double gamma_calc(double* u, const double &time_squared);
    KOKKOS_FUNCTION
    static double dot(double* v, double* u, const double &time_squared);
    KOKKOS_FUNCTION
    static void dot(double (*v)[D], double (*T)[D][D], const double &time_squared,
              double (*result)[D]);
    KOKKOS_FUNCTION
    static double get_shvDD(double* pi_diag, const double &time_squared);
    KOKKOS_FUNCTION
    static double get_LRF(const double &lab, const double &gamma, const double &time_squared);

    KOKKOS_FUNCTION
    static void evaluate_time_derivatives( Cabana::AoSoA<CabanaParticle, DeviceType, VECTOR_LENGTH> &particles );
    
    std::string name = "Israel-Stewart";                    // name associated to EoM

    //==========================================================================
    Matrix<double,D,D> dpidtsub_fun( hydrodynamic_info<D> & hi )
    {
      Matrix<double,D,D> vsub;

      for (unsigned int i=0; i<D; i++)
      for (unsigned int j=0; j<D; j++)
      for (unsigned int k=0; k<D; k++)
        vsub(i,j) += ( hi.u(i)*hi.pimin(j,k) + hi.u(j)*hi.pimin(i,k) )*hi.du_dt(k);

      return vsub;
    }

    //==========================================================================
    double Bsub_fun( hydrodynamic_info<D> & hi )
    {
      // make sure this quantity is set
      hi.uu = hi.u*hi.u;

      if ( !this->settingsPtr->using_shear )
        return 0.0;
      else
      {
        // these quantities will all be zero without shear
        mini( hi.pimin, hi.shv );
        hi.piu         = rowp1(0,hi.shv)*hi.u;
        hi.piutot      = hi.piu+transpose(hi.piu);
        double bsub = 0.0;
        double pig  = hi.shv(0,0)/hi.gamma_squared;

        for (unsigned int i=0; i<D; i++)
        for (unsigned int j=0; j<D; j++)
          bsub += hi.gradU(i,j) * ( hi.pimin(i,j) + pig*hi.uu(j,i)
                                - ( hi.piu(i,j) + hi.piu(j,i) ) / hi.gamma );
        return bsub;
      }
    }

    //==========================================================================
    Matrix<double,D,D> Msub_fun( hydrodynamic_info<D> & hi )
    {
      hi.piu  = rowp1(0,hi.shv)*hi.u;
      return hi.Agam2*hi.uu + hi.Ctot*hi.gamma*Imat - (1+4/3./hi.gamma_squared)*hi.piu
              + hi.dwdsT1*transpose(hi.piu) + hi.gamma*hi.pimin;
    }

    
    //==========================================================================
    
    
};
}
#endif