#ifndef EOM_DEFAULT_H
#define EOM_DEFAULT_H

#include "eom.h"
#include "hydrodynamic_info.h"
#include "matrix.h"
#include "settings.h"
#include "vector.h"

class EoM_default: public EquationsOfMotion
{

  public:
    EoM_default(){}
    virtual ~EoM_default(){}

    std::string name = "";                    // name associated to EoM

    // require all of these to be defined
    void compute_du_dt(){}
    void compute_dshv_dt(){}
    void compute_detasigma_dt(){}
    void compute_dBulk_dt(){}

    
    void evaluate_time_derivatives( hydrodynamic_info & hi )
    {
//cout << "t=: In " << __FILE__ << "::" << __LINE__ << endl;
      double gamt = 0.0, pre = 0.0, p1 = 0.0;
      if ( settingsPtr->using_shear )
      {
        gamt = 1.0/hi.gamma/hi.stauRelax;
        pre  = hi.eta_o_tau/hi.gamma;
        p1   = gamt - 4.0/3.0/hi.sigma*hi.dsigma_dt + 1.0/hi.t/3.0;
      }

      Vector<double,2> minshv   = rowp1(0, hi.shv);
      Matrix <double,2,2> partU = hi.gradU + transpose( hi.gradU );

      // set the Mass and the Force
      Matrix <double,2,2> M = hi.Msub;
      Vector<double,2> F    = hi.Btot*hi.u + hi.gradshear
                              - ( hi.gradP + hi.gradBulk + hi.divshear );

      //===============
      // shear contribution
      if ( settingsPtr->using_shear )
        F += pre*hi.v*partU + p1*minshv;

      double det = deter(M);

      Matrix <double,2,2> MI;
      MI(0,0) =  M(1,1)/det;
      MI(0,1) = -M(0,1)/det;
      MI(1,0) = -M(1,0)/det;
      MI(1,1) =  M(0,0)/det;

      //===============
      // compute acceleration
      hi.du_dt(0) = F(0) * MI(0,0) + F(1) * MI(0,1);
      hi.du_dt(1) = F(0) * MI(1,0) + F(1) * MI(1,1);

      //===============
      // define auxiliary variables
      double vduk               = inner( hi.v, hi.du_dt );
      Matrix <double,2,2> ulpi  = hi.u*colp1(0, hi.shv);
      Matrix <double,2,2> Ipi   = - 2.0*hi.eta_o_tau/3.0 * ( hi.Imat + hi.uu )
                                  + 4./3.*hi.pimin;

      //===============
      // "coordinate" divergence
      hi.div_u                   = (1./ hi.gamma)*inner( hi.u, hi.du_dt)
                                    - ( hi.gamma/ hi.sigma ) * hi.dsigma_dt;
      //===============
      // "covariant" divergence
      hi.bigtheta                = hi.div_u*hi.t+hi.gamma;

      //===============
      Matrix <double,2,2> sub   = hi.pimin + (hi.shv(0,0)/hi.g2)*hi.uu
                                  - 1./hi.gamma*hi.piutot;

      //===============
      if ( settingsPtr->using_shear )
        hi.inside                  = hi.t*( inner( -minshv+hi.shv(0,0)*hi.v, hi.du_dt )
                                      - con2(sub, hi.gradU) - hi.gamma*hi.t*hi.shv33 );

      hi.detasigma_dt            = 1./hi.sigma/hi.T*( -hi.bigPI*hi.bigtheta + hi.inside );


      // N.B. - ADD EXTRA TERMS FOR BULK EQUATION
      hi.dBulk_dt = ( -hi.zeta/hi.sigma*hi.bigtheta - hi.Bulk/hi.gamma )/hi.tauRelax;

      Matrix <double,2,2> ududt = hi.u*hi.du_dt;

      // N.B. - ADD READABLE TERM NAMES
      if ( settingsPtr->using_shear )
        hi.dshv_dt                 = - gamt*( hi.pimin + hi.setas*partU )
                                 - hi.eta_o_tau*( ududt + transpose(ududt) )
                                 + hi.dpidtsub + hi.sigl*Ipi
                                 - vduk*( ulpi + transpose(ulpi) + (1/hi.gamma)*Ipi );
std::cout << "t=CHECK: " << hi.ID << "   "
              << hi.t << "   "
              << hi.dBulk_dt << "   "
              << hi.detasigma_dt << "   "
              << hi.du_dt << "   "
              << hi.dshv_dt << "\n";

if (hi.t > 1.5) exit(8);

    }

};


#endif