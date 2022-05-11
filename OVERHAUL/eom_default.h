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


    //==========================================================================
    Matrix<double,2,2> dpidtsub_fun( hydrodynamic_info & hi )
    {
      Matrix<double,2,2> vsub;

      for (int i=0; i<=1; i++)
      for (int j=0; j<=1; j++)
      for (int k=0; k<=1; k++)
        vsub(i,j) += ( hi.u(i)*hi.pimin(j,k) + hi.u(j)*hi.pimin(i,k) )*hi.du_dt(k);

      return vsub;
    }

    //==========================================================================
    double Bsub_fun()
    {
      // make sure this quantity is set
      uu = u*u;

      if ( !settingsPtr->using_shear )
        return 0.0;
      else
      {
        // these quantities will all be zero without shear
        mini( pimin, shv );
        piu         = rowp1(0,shv)*u;
        piutot      = piu+transpose(piu);
        double bsub = 0.0;
        double pig  = shv(0,0)/g2;

        for (int i=0; i<=1; i++)
        for (int j=0; j<=1; j++)
          bsub += gradU(i,j) * ( pimin(i,j) + pig*uu(j,i)
                                - ( piu(i,j) + piu(j,i) ) / gamma );
        return bsub;
      }
    }

    //==========================================================================
    Matrix<double,2,2> Msub_fun()
    {
      piu  = rowp1(0,shv)*u;
      return Agam2*uu + Ctot*gamma*Imat - (1+4/3./g2)*piu
              + dwdsT1*transpose(piu) + gamma*pimin;
    }

    
    //==========================================================================
    void evaluate_time_derivatives( hydrodynamic_info & hi )
    {
      // PREVIOUSLY DONE IN UPDATE_DSIGMA_DT
      hi.dsigma_dt = -hi.sigma * ( hi.gradV(0,0) + hi.gradV(1,1) );






      // PREVIOUSLY DONE IN UPDATE_FLUID_VARIABLES
      hi.g2           = hi.gamma*hi.gamma;
      hi.g3           = hi.gamma*hi.g2;
      hi.gt           = hi.gamma*hi.t;
      double dwdsT    = hi.dwds/hi.T;
      hi.dwdsT1       = 1 - hi.dwds/hi.T;
      hi.sigl         = hi.dsigma_dt/hi.sigma - 1/hi.t;
      hi.gradU        = hi.gamma*hi.gradV+hi.g3*(hi.v*(hi.v*hi.gradV));
      hi.bigPI        = hi.Bulk*hi.sigma/hi.gt;
      hi.C            = hi.w+ hi.bigPI;

      hi.eta_o_tau    = (settingsPtr->using_shear) ? hi.setas/hi.stauRelax : 0.0;

      hi.Agam         = hi.w - hi.dwds*( hi.s + hi.bigPI/hi.T ) - hi.zeta/hi.tauRelax
                        - hi.dwdB*hi.rhoB - hi.dwdS*hi.rhoS - hi.dwdQ*hi.rhoQ;

      hi.Agam2        = ( hi.Agam - hi.eta_o_tau/3.0 - hi.dwdsT1*hi.shv(0,0) ) / hi.gamma;
      hi.Ctot         = hi.C + hi.eta_o_tau*(1.0/hi.g2-1.0);


      hi.Btot         = ( hi.Agam*hi.gamma + 2.0*hi.eta_o_tau/3.0*hi.gamma )*hi.sigl
                          + hi.bigPI/hi.tauRelax
                          + hi.dwdsT*( hi.gt*hi.shv33 + Bsub_fun() );
      hi.check        = hi.sigl;









      // THIS IS THE ORIGINAL PART IN TIME DERIVATIVES
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

      hi.dpidtsub = dpidtsub_fun(hi);

      // N.B. - ADD READABLE TERM NAMES
      if ( settingsPtr->using_shear )
        hi.dshv_dt                 = - gamt*( hi.pimin + hi.setas*partU )
                                 - hi.eta_o_tau*( ududt + transpose(ududt) )
                                 + hi.dpidtsub + hi.sigl*Ipi
                                 - vduk*( ulpi + transpose(ulpi) + (1/hi.gamma)*Ipi );
//std::cout << "t=CHECK: " << hi.ID << "   "
//              << hi.t << "   "
//              << hi.dBulk_dt << "   "
//              << hi.detasigma_dt << "   "
//              << hi.du_dt << "   "
//              << hi.dshv_dt << "\n";
//
//std::cout << "CHECK dshv_dt: " << hi.ID << "   "
//              << hi.t << "   "
//              << hi.dshv_dt << "   "
//              << - gamt*( hi.pimin + hi.setas*partU ) << "   "
//              << - hi.eta_o_tau*( ududt + transpose(ududt) ) << "   "
//              << hi.dpidtsub + hi.sigl*Ipi << "   "
//              << - vduk*( ulpi + transpose(ulpi) + (1/hi.gamma)*Ipi ) << "\n";
//
//if (hi.t > 1.1) exit(8);

    }

};


#endif