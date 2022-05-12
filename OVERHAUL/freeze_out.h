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

    // freeze out struct
    struct FRZ
    {
      double t = 0.0, s = 0.0, e = 0.0, rhoB = 0.0, rhoS = 0.0, rhoQ = 0.0,
             T = 0.0, muB = 0.0, muS = 0.0, muQ = 0.0, theta = 0.0, bulk = 0.0,
             sigma = 0.0, shear33 = 0.0, inside = 0.0;
      Vector<double,2> r, u, gradP;
      Matrix<double,3,3> shear;
    };

    vector<FRZ> frz1;
    vector<FRZ> frz2;
    vector<FRZ> fback;
    vector<FRZ> fback2;
    vector<FRZ> fback3;
    vector<FRZ> fback4;


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
      if ( freeze_out_mode != "constant_energy_density" )
      {
        cerr << "This freeze out mode is not supported.  Please try a different one.\n";
        exit(8);
      }

      if ( systemPtr == nullptr )
      {
        cerr << "You need to provide the location of the SystemState object.\n";
        exit(8);
      }

      // resize vectors to contain all particles
      frz1.resize( systemPtr->particles.size() );
      frz2.resize( systemPtr->particles.size() );
      fback.resize( systemPtr->particles.size() );
      fback2.resize( systemPtr->particles.size() );
      fback3.resize( systemPtr->particles.size() );
      fback4.resize( systemPtr->particles.size() );
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
          auto & p_frz2  = frz2[p.ID];
          p_frz2.r       = p.r;
          p_frz2.u       = p.hydro.u;
          p_frz2.sigma   = p.hydro.sigma;
          p_frz2.T       = p.T();
          p_frz2.muB     = p.muB();
          p_frz2.muS     = p.muS();
          p_frz2.muQ     = p.muQ();
          p_frz2.e       = p.e();
          p_frz2.rhoB    = p.rhoB();
          p_frz2.rhoS    = p.rhoS();
          p_frz2.rhoQ    = p.rhoQ();
          p_frz2.bulk    = p.hydro.bigPI;
          p_frz2.theta   = p.hydro.div_u + p.hydro.gamma/systemPtr->t;
          p_frz2.gradP   = p.hydro.gradP;
          p_frz2.shear   = p.hydro.shv;
          p_frz2.shear33 = p.hydro.shv33;
          p_frz2.inside  = p.hydro.inside;
        }

      }
      else if (systemPtr->frzc==1)
      {
        systemPtr->taup = systemPtr->t;
        systemPtr->frzc = 2;
        for (auto & p : systemPtr->particles)
        {
          auto & p_frz1  = frz1[p.ID];
          p_frz1.r       = p.r;
          p_frz1.u       = p.hydro.u;
          p_frz1.sigma   = p.hydro.sigma;
          p_frz1.T       = p.T();
          p_frz1.muB     = p.muB();
          p_frz1.muS     = p.muS();
          p_frz1.muQ     = p.muQ();
          p_frz1.e       = p.e();
          p_frz1.rhoB    = p.rhoB();
          p_frz1.rhoS    = p.rhoS();
          p_frz1.rhoQ    = p.rhoQ();
          p_frz1.bulk    = p.hydro.bigPI;
          p_frz1.theta   = p.hydro.div_u + p.hydro.gamma/systemPtr->t;
          p_frz1.gradP   = p.hydro.gradP;
          p_frz1.shear   = p.hydro.shv;
          p_frz1.shear33 = p.hydro.shv33;
          p_frz1.inside  = p.hydro.inside;
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
              int p_ID = p.ID;
              fback4[p_ID] = fback2[p_ID];
              fback3[p_ID] = fback[p_ID];
              fback2[p_ID] = frz2[p_ID];
              fback[p_ID]  = frz1[p_ID];
            }
            else if ( p.btrack == 0 )
            {
              if ( p.fback.gradP(0) != 0 )
              {
                int p_ID = p.ID;
                frz2[p_ID] = fback2[p_ID];
                frz1[p_ID] = fback[p_ID];
              }
              else
              {
                int p_ID = p.ID;
                frz2[p_ID] = fback4[p_ID];
                frz1[p_ID] = fback3[p_ID];
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
          auto & p_frz1  = frz1[p.ID];
          auto & p_frz2  = frz2[p.ID];
          p_frz2         = p_frz1;

          p_frz1.r       = p.r;
          p_frz1.u       = p.hydro.u;
          p_frz1.sigma   = p.hydro.sigma;
          p_frz1.T       = p.T();
          p_frz1.muB     = p.muB();
          p_frz1.muS     = p.muS();
          p_frz1.muQ     = p.muQ();
          p_frz1.e       = p.e();
          p_frz1.rhoB    = p.rhoB();
          p_frz1.rhoS    = p.rhoS();
          p_frz1.rhoQ    = p.rhoQ();
          p_frz1.bulk    = p.hydro.bigPI ;
          p_frz1.theta   = p.hydro.div_u+p.hydro.gamma/systemPtr->t;
          p_frz1.gradP   = p.hydro.gradP;
          p_frz1.shear   = p.hydro.shv;
          p_frz1.shear33 = p.hydro.shv33;
          p_frz1.inside  = p.hydro.inside;
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
        if ( abs( frz1[i].e - systemPtr->efcheck ) < abs( frz2[i].e - systemPtr->efcheck ) )
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

          systemPtr->rsub[j]       = frz1[i].r;
          systemPtr->uout[j]       = frz1[i].u;
          systemPtr->bulksub[j]    = frz1[i].bulk;
          systemPtr->shearsub[j]   = frz1[i].shear;
          systemPtr->shear33sub[j] = frz1[i].shear33;

          gradPsub      = frz1[i].gradP;
          inside        = frz1[i].inside;
          sigsub        = frz1[i].sigma;
          thetasub      = frz1[i].theta;
          systemPtr->Tfluc[j]      = frz1[i].T;             // replace with e
        }
        else if ( swit == 2 )
        {
          if ( p.btrack != -1 )
            systemPtr->tlist[j]    = systemPtr->taupp;
          else
            systemPtr->tlist[j]    = systemPtr->taupp - systemPtr->dt;

          systemPtr->rsub[j]       = frz2[i].r;
          systemPtr->uout[j]       = frz2[i].u;
          systemPtr->bulksub[j]    = frz2[i].bulk;
          systemPtr->shearsub[j]   = frz2[i].shear;
          systemPtr->shear33sub[j] = frz2[i].shear33;

          gradPsub      = frz2[i].gradP;
          inside        = frz2[i].inside;
          sigsub        = frz2[i].sigma;
          thetasub      = frz2[i].theta;
          systemPtr->Tfluc[j]      = frz2[i].T;           // replace with e
        }
        else
        {
          cout << __PRETTY_FUNCTION__ << ": Not at freeze-out temperature" << endl;
        }

        systemPtr->sFO[j]       = eosPtr->s_terms_T( systemPtr->Tfluc[j] );  // replace with e, BSQ
        systemPtr->gsub[j]      = sqrt( Norm2(systemPtr->uout[j]) + 1 );


        sigsub      /= systemPtr->gsub[j]*systemPtr->tlist[j];
        systemPtr->swsub[j]     = p.norm_spec.s/sigsub;

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
          cout << frz1[i].gradP << " " << frz2[i].gradP << endl;
          cout << frz1[i].T*197.3<< " " << frz2[i].T*197.3 << endl;
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
          cout << frz1[i].T*0.1973<< " " << frz2[i].T*0.1973<< endl;
        }

        systemPtr->sFO[j]   *= pow(systemPtr->Tfluc[j]*0.1973, 3);
        systemPtr->Tfluc[j] *= 0.1973;

      }

      systemPtr->cf = curfrz;
    }


    //==============================================================================
    void frzcheck( Particle & p, double tin, int &count, int N )
    {
      int p_ID = p.ID;
      if ( p.Freeze == 0 )
      {
        if ( p.e() <= p.efcheck )
        {
          p.Freeze = 1;
          frz2[p_ID].t = tin;
        }
      }
      else if ( p.Freeze == 1 )
      {
        if ( p.btrack == -1 )
        {
          count += 1;
          p.Freeze = 3;
          frz1[p_ID].t = tin;
        }
        else if ( p.e()>frz1[p_ID].e )
        {
          p.Freeze = 1;
          frz2[p_ID].t = tin;
        }
        else if( p.e() <= p.efcheck )
        {
          count += 1;
          p.Freeze = 3;
          frz1[p_ID].t = tin;
        }
        else
        {
          p.Freeze=0;
        }
      }

    }


};

typedef std::shared_ptr<FreezeOut> pFreezeOut;  // smart pointer to freeze out object

#endif