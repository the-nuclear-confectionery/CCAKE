#ifndef FREEZEOUT_H
#define FREEZEOUT_H

#include <algorithm>
#include <string>
#include <vector>

#include "eos.h"
#include "particle.h"
#include "system_state.h"

using std::string;

class FreezeOut
{
  private:
    
    EquationOfState * eosPtr  = nullptr;
    Settings * settingsPtr    = nullptr;
    SystemState * systemPtr   = nullptr;

    string freeze_out_mode = "constant_energy_density";

    double freeze_out_threshold = 0.0;  // units depend on what mode is used

  public:

    void set_EquationOfStatePtr(EquationOfState * eosPtr_in) { eosPtr = eosPtr_in; }
    void set_SettingsPtr(Settings * settingsPtr_in) { settingsPtr = settingsPtr_in; }
    void set_SystemStatePtr(SystemState * systemPtr_in) { systemPtr = systemPtr_in; }

    // set things up
    void initialize( double fo_threshold_in,
                     const string & fo_mode = "constant_energy_density" )
    {
      freeze_out_threshold = fo_threshold_in;
      freeze_out_mode      = fo_mode;
    }

    ////////////////////////////////////////////////////////////////////////////
    void bsqsvfreezeout(int curfrz)
    {
      cout << "CHECK BSQSVFREEZEOUT: " << systemPtr->frzc
            << "   " << systemPtr->tau << "   " << systemPtr->taup
            << "   " << systemPtr->taupp << "   " << systemPtr->cfon << endl;

      if (systemPtr->frzc==0)
      {
        systemPtr->taupp = systemPtr->t;
        systemPtr->frzc  = 1;
        for (auto & p : systemPtr->particles)
        {
          p.frz2.r       = p.r;
          p.frz2.u       = p.hydro.u;
          p.frz2.sigma   = p.hydro.sigma;
          p.frz2.T       = p.T();
          p.frz2.muB     = p.muB();
          p.frz2.muS     = p.muS();
          p.frz2.muQ     = p.muQ();
          p.frz2.e       = p.e();
          p.frz2.rhoB    = p.rhoB();
          p.frz2.rhoS    = p.rhoS();
          p.frz2.rhoQ    = p.rhoQ();
          p.frz2.bulk    = p.hydro.bigPI;
          p.frz2.theta   = p.hydro.div_u + p.hydro.gamma/systemPtr->t;
          p.frz2.gradP   = p.hydro.gradP;
          p.frz2.shear   = p.hydro.shv;
          p.frz2.shear33 = p.hydro.shv33;
          p.frz2.inside  = p.hydro.inside;
        }

      }
      else if (systemPtr->frzc==1)
      {
        systemPtr->taup = systemPtr->t;
        systemPtr->frzc = 2;
        for (auto & p : systemPtr->particles)
        {
          p.frz1.r       = p.r;
          p.frz1.u       = p.hydro.u;
          p.frz1.sigma   = p.hydro.sigma;
          p.frz1.T       = p.T();
          p.frz1.muB     = p.muB();
          p.frz1.muS     = p.muS();
          p.frz1.muQ     = p.muQ();
          p.frz1.e       = p.e();
          p.frz1.rhoB    = p.rhoB();
          p.frz1.rhoS    = p.rhoS();
          p.frz1.rhoQ    = p.rhoQ();
          p.frz1.bulk    = p.hydro.bigPI;
          p.frz1.theta   = p.hydro.div_u + p.hydro.gamma/systemPtr->t;
          p.frz1.gradP   = p.hydro.gradP;
          p.frz1.shear   = p.hydro.shv;
          p.frz1.shear33 = p.hydro.shv33;
          p.frz1.inside  = p.hydro.inside;
        }

        systemPtr->divTtemp.resize( curfrz );
        systemPtr->divT.resize( curfrz );
        systemPtr->gsub.resize( curfrz );
        systemPtr->uout.resize( curfrz );
        systemPtr->swsub.resize( curfrz );
        systemPtr->bulksub.resize( curfrz );
        systemPtr->shearsub.resize( curfrz );
        systemPtr->shear33sub.resize( curfrz );
        systemPtr->tlist.resize( curfrz );
        systemPtr->rsub.resize( curfrz );

        if ( curfrz > 0 )
          bsqsvinterpolate( curfrz );
        else
          systemPtr->cf = 0;
      }
      else
      {
        int i_local = 0;
        for (auto & p : systemPtr->particles)
        {
          if ( p.Freeze < 4 )
          {
            if ( ( p.btrack <= 3 ) && ( p.btrack > 0 ) )
            {
              p.fback4 = p.fback2;
              p.fback3 = p.fback;
              p.fback2 = p.frz2;
              p.fback  = p.frz1;
            }
            else if ( p.btrack == 0 )
            {
              if ( p.fback.gradP(0) != 0 )
              {
                p.frz2 = p.fback2;
                p.frz1 = p.fback;
              }
              else
              {
                p.frz2 = p.fback4;
                p.frz1 = p.fback3;
                cout << "back second"  << endl;
              }


              curfrz++;
              systemPtr->list.push_back( i_local );
              p.Freeze = 4;
              p.btrack = -1;
            }
          }

          i_local++;
        }

        systemPtr->tau = systemPtr->t;

        // resize vectors
        systemPtr->divTtemp.resize( curfrz );
        systemPtr->divT.resize( curfrz );
        systemPtr->gsub.resize( curfrz );
        systemPtr->uout.resize( curfrz );
        systemPtr->swsub.resize( curfrz );
        systemPtr->bulksub.resize( curfrz );
        systemPtr->shearsub.resize( curfrz );
        systemPtr->shear33sub.resize( curfrz );
        systemPtr->tlist.resize( curfrz );
        systemPtr->rsub.resize( curfrz );

        if ( curfrz > 0 )
          bsqsvinterpolate( curfrz );
        else
          systemPtr->cf = 0;


        //sets up the variables for the next time step
        for (auto & p : systemPtr->particles)
        {
          p.frz2         = p.frz1;

          p.frz1.r       = p.r;
          p.frz1.u       = p.hydro.u;
          p.frz1.sigma   = p.hydro.sigma;
          p.frz1.T       = p.T();
          p.frz1.muB     = p.muB();
          p.frz1.muS     = p.muS();
          p.frz1.muQ     = p.muQ();
          p.frz1.e       = p.e();
          p.frz1.rhoB    = p.rhoB();
          p.frz1.rhoS    = p.rhoS();
          p.frz1.rhoQ    = p.rhoQ();
          p.frz1.bulk    = p.hydro.bigPI ;
          p.frz1.theta   = p.hydro.div_u+p.hydro.gamma/systemPtr->t;
          p.frz1.gradP   = p.hydro.gradP;
          p.frz1.shear   = p.hydro.shv;
          p.frz1.shear33 = p.hydro.shv33;
          p.frz1.inside  = p.hydro.inside;
        }

        systemPtr->taupp = systemPtr->taup;
        systemPtr->taup  = systemPtr->tau;
      }

      systemPtr->cfon = 0;
    }


    //////////////////////////////////////////////////////////////////////////////
    void bsqsvinterpolate(int curfrz)
    {
      systemPtr->sFO.resize( curfrz, 0 );
      systemPtr->Tfluc.resize( curfrz, 0 );

      for (int j=0; j<curfrz; j++)
      {
        int i    = systemPtr->list[j];
        auto & p = systemPtr->particles[i];

        int swit = 0;
        if ( abs( p.frz1.e - systemPtr->efcheck ) < abs( p.frz2.e - systemPtr->efcheck ) )
          swit   = 1;
        else
          swit   = 2;

        double sigsub = 0.0, thetasub = 0.0, inside = 0.0;
        Vector<double,2> gradPsub;
        if ( swit == 1 )
        {
          if ( p.btrack != -1 )
            systemPtr->tlist[j]    = systemPtr->taup;
          else
            systemPtr->tlist[j]    = systemPtr->taup - systemPtr->dt;

          systemPtr->rsub[j]       = p.frz1.r;
          systemPtr->uout[j]       = p.frz1.u;
          systemPtr->bulksub[j]    = p.frz1.bulk;
          systemPtr->shearsub[j]   = p.frz1.shear;
          systemPtr->shear33sub[j] = p.frz1.shear33;

          gradPsub      = p.frz1.gradP;
          inside        = p.frz1.inside;
          sigsub        = p.frz1.sigma;
          thetasub      = p.frz1.theta;
          systemPtr->Tfluc[j]      = p.frz1.T;             // replace with e
        }
        else if ( swit == 2 )
        {
          if ( p.btrack != -1 )
            systemPtr->tlist[j]    = systemPtr->taupp;
          else
            systemPtr->tlist[j]    = systemPtr->taupp - systemPtr->dt;

          systemPtr->rsub[j]       = p.frz2.r;
          systemPtr->uout[j]       = p.frz2.u;
          systemPtr->bulksub[j]    = p.frz2.bulk;
          systemPtr->shearsub[j]   = p.frz2.shear;
          systemPtr->shear33sub[j] = p.frz2.shear33;

          gradPsub      = p.frz2.gradP;
          inside        = p.frz2.inside;
          sigsub        = p.frz2.sigma;
          thetasub      = p.frz2.theta;
          systemPtr->Tfluc[j]      = p.frz2.T;           // replace with e
        }
        else
        {
          cout << __PRETTY_FUNCTION__ << ": Not at freeze-out temperature" << endl;
        }

        systemPtr->sFO[j]       = eosPtr->s_terms_T( systemPtr->Tfluc[j] );  // replace with e, BSQ
        systemPtr->gsub[j]      = sqrt( Norm2(systemPtr->uout[j]) + 1 );


        sigsub      /= systemPtr->gsub[j]*systemPtr->tlist[j];
        systemPtr->swsub[j]     = p.sigmaweight/sigsub;

        systemPtr->divT[j]      = (1.0/systemPtr->sFO[j])*gradPsub;
        systemPtr->divTtemp[j]  = -(1.0/(systemPtr->gsub[j]*systemPtr->sFO[j]))
                          *( systemPtr->cs2 * (systemPtr->wfz+systemPtr->bulksub[j]) * thetasub
                            - systemPtr->cs2*inside+inner(systemPtr->uout[j], gradPsub) );
    //THIS NEEDS TO BE RESET


        double insub = systemPtr->divTtemp[j]*systemPtr->divTtemp[j] - Norm2(systemPtr->divT[j]);
        double norm  = -sqrt(abs(insub));
        systemPtr->divTtemp[j] /= norm;
        systemPtr->divT[j]      = (1.0/norm)*systemPtr->divT[j];


        if ( systemPtr->divTtemp[j] == 1 )
        {
          cout << "track sph=" << p.btrack << " " << i << endl;
          cout << systemPtr->divTtemp[j] << " " << systemPtr->divT[j] << " " << norm << endl;
          cout << gradPsub << " " << thetasub << endl;
          cout << systemPtr->tlist[j] << " " << p.r << endl;
          cout << p.frz1.gradP << " " << p.frz2.gradP << endl;
          cout << p.frz1.T*197.3<< " " << p.frz2.T*197.3 << endl;
          getchar();
        }

        systemPtr->avgetasig += systemPtr->sFO[j]/sigsub;

        if(isnan(systemPtr->divTtemp[j]))
        {
          cout << "divtemp" << endl;
          cout << systemPtr->divTtemp[j] << " " << systemPtr->divT[j] << " " << norm << endl;
          cout << gradPsub << " " << thetasub << endl;
          cout << systemPtr->bulksub[j] << endl;
          cout << systemPtr->gsub[j] << endl;
          cout << systemPtr->tlist[j] << " " << p.r << endl;
          cout << p.frz1.T*0.1973<< " " << p.frz2.T*0.1973<< endl;
        }

        systemPtr->sFO[j]   *= pow(systemPtr->Tfluc[j]*0.1973, 3);
        systemPtr->Tfluc[j] *= 0.1973;

      }

      systemPtr->cf = curfrz;
    }
};

typedef std::shared_ptr<FreezeOut> pFreezeOut;  // smart pointer to freeze out object

#endif