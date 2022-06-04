#ifndef FREEZEOUT_H
#define FREEZEOUT_H

#include <algorithm>
#include <string>
#include <vector>

#include "eos.h"
#include "input_output.h"
#include "particle.h"
#include "system_state.h"

using std::string;

//==============================================================================
// A DICTIONARY FOR CONFUSINGLY NAMED QUANTITIES
//------------------------------------------------------------------------------
// tau    - current tau
// taup   - tau of previous timestep
// taupp  - tau of two timesteps ago
// Freeze - current freeze out status of particle (0 == freeze-out not begun,
//                                                 1 == freeze-out begun,
//                                                 3 == freeze-out basically done,
//                                                 4 == completely frozen out)
// frz1   - snapshot of particle state at taup
// frz2   - snapshot of particle state at taupp
// fback1 - snapshot of particle state at taup
// fback2 - snapshot of particle state at taupp
// fback3 - snapshot of particle state at tauppp(?)
// fback4 - snapshot of particle state at taupppp(?)
// cf     - basically the same as curfrz
// btrack - essentially counts the number of nearest neighbors a particle has
//==============================================================================

class FreezeOut
{
  friend class InputOutput;

  private:
    
//    EquationOfState * eosPtr  = nullptr;
    Settings * settingsPtr    = nullptr;
    SystemState * systemPtr   = nullptr;

    string freeze_out_mode = "constant_energy_density";

    double freeze_out_threshold = 0.0;  // units depend on what mode is used

    // freeze out struct
    struct FRZ
    {
      string eos_name = "";
      double t = 0.0, s = 0.0, e = 0.0, rhoB = 0.0, rhoS = 0.0, rhoQ = 0.0,
             T = 0.0, muB = 0.0, muS = 0.0, muQ = 0.0, theta = 0.0, bulk = 0.0,
             sigma = 0.0, shear33 = 0.0, inside = 0.0;
      Vector<double,2> r, u, gradP;
      Matrix<double,3,3> shear;
    };


    //==========================================================================
    // MOVE TO THEIR OWN CLASS???
    // freeze-out related quantities
    int cf            = 0;
    int frzc          = 0;
    double avgetasig  = 0;
    double cs2        = 0;
    double tau        = 0;
    double taup       = 0;
    double taupp      = 0;
    double wfz        = 0;

    vector<double> divTtemp;
    vector<double> gsub;
    vector<double> bulksub;
    vector<double> swsub;
    vector<double> shear33sub;
    vector<double> tlist;
    vector<double> sFO;         //entropy at freezeout
    vector<double> Tfluc, muBfluc, muSfluc, muQfluc;

    vector<string> eosname;

    vector< Vector<double,2> > divT;
    vector< Vector<double,2> > rsub;
    vector< Vector<double,2> > uout;
    vector< Matrix<double,3,3> > shearsub;
    //==========================================================================    


    vector<FRZ> frz1;
//    vector<FRZ> frz2;
    vector<FRZ> fback;
    vector<FRZ> fback2;
    vector<FRZ> fback3;
    vector<FRZ> fback4;

  public:

    // for some dumb reason, this needs to be public
    // (must be accessed in system state)
    vector<FRZ> frz2;

//    void set_EquationOfStatePtr(EquationOfState * eosPtr_in) { eosPtr = eosPtr_in; }
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








    //==============================================================================
    void check_freeze_out_status( Particle & p, double tin, int &count, int N )
    {
      int p_ID = p.ID;
      if ( p.Freeze == 0 )              // if particle has not yet started freezing out
      {
        if ( p.e() <= p.efcheck )       //  * if it's energy is below eFO
        {
          p.Freeze = 1;                 //    - start freeze-out
          frz2[p_ID].t = tin;           //    - note the time as the most distant
        }
      }
      else if ( p.Freeze == 1 )         // otherwise, if freeze-out has begun
      {
        if ( p.btrack == -1 )           //  * if particle has no neighbors
        {
          count += 1;                   //    - increment currently freezing out count
          p.Freeze = 3;                 //    - set Freeze = 3
          frz1[p_ID].t = tin;           //    - note the time as previous
        }
        else if ( p.e()>frz1[p_ID].e )  //  * otherwise, if the particle's energy
        {                               //    has increased since the previous timestep
          p.Freeze = 1;                 //    - keep Freeze = 1
          frz2[p_ID].t = tin;           //    - set penultimate time
        }
        else if( p.e() <= p.efcheck )   //  * otherwise, if e is below eFO
        {
          count += 1;                   //    - increment currently freezing out count
          p.Freeze = 3;                 //    - set Freeze = 3
          frz1[p_ID].t = tin;           //    - note the time as previous
        }
        else                            //  * otherwise,
        {
          p.Freeze=0;                   //    - particle is not ready to freeze-out
        }
      }

    }









    ////////////////////////////////////////////////////////////////////////////
    void bsqsvfreezeout(int curfrz)
    {
      if (frzc==0) // true in first timestep only
      {
        taupp = systemPtr->t;
        frzc  = 1;
cout << "Check sizes: " << frz2.size() << "   " << systemPtr->particles.size() << endl;
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
          p_frz2.s       = p.s();
          p_frz2.rhoB    = p.rhoB();
          p_frz2.rhoS    = p.rhoS();
          p_frz2.rhoQ    = p.rhoQ();
          p_frz2.bulk    = p.hydro.bigPI;
          p_frz2.theta   = p.hydro.div_u + p.hydro.gamma/systemPtr->t;
          p_frz2.gradP   = p.hydro.gradP;
          p_frz2.shear   = p.hydro.shv;
          p_frz2.shear33 = p.hydro.shv33;
          p_frz2.inside  = p.hydro.inside;
//cout << __FUNCTION__ << ":" << __LINE__ << ": " << endl;
//cout << "check 1" << endl;
//cout << p_frz2.eos_name << endl;
//cout << "check 2" << endl;
//cout << p.get_current_eos_name() << endl;
//cout << "checks passed" << endl;
//cout << __FUNCTION__ << ":" << __LINE__ << ": " << endl;
          p_frz2.eos_name  = p.get_current_eos_name();
        }

      }
      else if (frzc==1) // true in second timestep only
      {
        taup = systemPtr->t;
        frzc = 2;
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
          p_frz1.s       = p.s();
          p_frz1.rhoB    = p.rhoB();
          p_frz1.rhoS    = p.rhoS();
          p_frz1.rhoQ    = p.rhoQ();
          p_frz1.bulk    = p.hydro.bigPI;
          p_frz1.theta   = p.hydro.div_u + p.hydro.gamma/systemPtr->t;
          p_frz1.gradP   = p.hydro.gradP;
          p_frz1.shear   = p.hydro.shv;
          p_frz1.shear33 = p.hydro.shv33;
          p_frz1.inside  = p.hydro.inside;
          p_frz1.eos_name  = p.get_current_eos_name();
        }

        divTtemp.resize( curfrz );
        divT.resize( curfrz );
        gsub.resize( curfrz );
        uout.resize( curfrz );
        swsub.resize( curfrz );
        bulksub.resize( curfrz );
        shearsub.resize( curfrz );
        shear33sub.resize( curfrz );
        tlist.resize( curfrz );
        rsub.resize( curfrz );

        if ( curfrz > 0 )
          bsqsvinterpolate( curfrz );
        else
          cf = 0;
      }
      else // true in third timestep and afterward
      {
        int i_local = 0;
        for (auto & p : systemPtr->particles)
        {
          int p_ID = p.ID;
          if ( p.Freeze < 4 )                             // particle not fully frozen out
          {
            if ( ( p.btrack <= 3 ) && ( p.btrack > 0 ) )  // particle has 1-3 nearest neighbors
            {
              fback4[p_ID] = fback2[p_ID];                // shift particle snapshots
              fback3[p_ID] = fback[p_ID];                 // at two past timesteps
              fback2[p_ID] = frz2[p_ID];                  // back another two timesteps
              fback[p_ID]  = frz1[p_ID];                  // (thus saving four previous
                                                          //  timesteps in total)
            }
            else if ( p.btrack == 0 )                     // particle has NO nearest neighbors
            {
              if ( fback[p_ID].gradP(0) != 0 )            // if particle still had neighbors
              {                                           // at the previous timestep,
                frz2[p_ID] = fback2[p_ID];                // REUSE last two timesteps
                frz1[p_ID] = fback[p_ID];                 // for freeze out interpolation
              }
              else
              {
                frz2[p_ID] = fback4[p_ID];                // if particle did NOT still
                frz1[p_ID] = fback3[p_ID];                // have neighbors at previous
                cout << "back second"  << endl;           // timestep, RECYCLE 3-4 timesteps ago
              }                                           // for freeze out interpolation

              //================================================================
              // !!!! N.B.: CURFRZ IS NOT ACTUALLY INCREMENTED HERE WITHOUT & ABOVE !!!!
              curfrz++;                                   // freeze-out the particle out
              systemPtr->list.push_back( i_local );       // list stores currently freezing
              p.Freeze = 4;                               // out particle indices
              p.btrack = -1;                              // don't let it "freeze in"?
            }
          }

          i_local++;
        }

        tau = systemPtr->t;

        // resize vectors
        divTtemp.resize( curfrz );
        divT.resize( curfrz );
        gsub.resize( curfrz );
        uout.resize( curfrz );
        swsub.resize( curfrz );
        bulksub.resize( curfrz );
        shearsub.resize( curfrz );
        shear33sub.resize( curfrz );
        tlist.resize( curfrz );
        rsub.resize( curfrz );

        if ( curfrz > 0 )
          bsqsvinterpolate( curfrz );
        else
          cf = 0;


        //sets up the variables for the next time step
        for (auto & p : systemPtr->particles)
        {
          auto & p_frz1  = frz1[p.ID];
          auto & p_frz2  = frz2[p.ID];
          p_frz2         = p_frz1;  // reset snapshots appropriately

          p_frz1.r       = p.r;
          p_frz1.u       = p.hydro.u;
          p_frz1.sigma   = p.hydro.sigma;
          p_frz1.T       = p.T();
          p_frz1.muB     = p.muB();
          p_frz1.muS     = p.muS();
          p_frz1.muQ     = p.muQ();
          p_frz1.e       = p.e();
          p_frz1.s       = p.s();
          p_frz1.rhoB    = p.rhoB();
          p_frz1.rhoS    = p.rhoS();
          p_frz1.rhoQ    = p.rhoQ();
          p_frz1.bulk    = p.hydro.bigPI ;
          p_frz1.theta   = p.hydro.div_u+p.hydro.gamma/systemPtr->t;
          p_frz1.gradP   = p.hydro.gradP;
          p_frz1.shear   = p.hydro.shv;
          p_frz1.shear33 = p.hydro.shv33;
          p_frz1.inside  = p.hydro.inside;
          p_frz1.eos_name  = p.get_current_eos_name();
        }

        taupp = taup;
        taup  = tau;
      }

      systemPtr->do_freeze_out = false;
    }







    //==========================================================================
    // THIS FUNCTION APPLIES ONLY TO PARTICLES WHICH ARE CURRENTLY IN THE PROCESS
    // OF FREEZING OUT
    void bsqsvinterpolate(int curfrz)
    {
      //cout << "Entered " << __FUNCTION__ << " with " << curfrz << " particles" << endl;
      sFO.resize( curfrz, 0 );
      Tfluc.resize( curfrz, 0 );
      muBfluc.resize( curfrz, 0 );
      muSfluc.resize( curfrz, 0 );
      muQfluc.resize( curfrz, 0 );
      eosname.resize( curfrz );

      // for all CURRENTLY FREEZING OUT PARTICLES
      for (int j=0; j<curfrz; j++)
      {
        int i    = systemPtr->list[j];
        auto & p = systemPtr->particles[i];

        // decide whether the particle was closer to freeze out at the previous
        // timestep or the one before that
        int swit = 0;
        if ( abs( frz1[i].e - systemPtr->efcheck ) < abs( frz2[i].e - systemPtr->efcheck ) )
          swit   = 1;
        else
          swit   = 2;

        double sigsub = 0.0, thetasub = 0.0, inside = 0.0;
        Vector<double,2> gradPsub;
        if ( swit == 1 )  // if particle was closer to freeze-out at last timestep
        {
          if ( p.btrack != -1 )                 // if particle had neighbors
            tlist[j]    = taup;                 // at previous timestep
          else
            tlist[j]    = taup - systemPtr->dt; // otherwise, go back two timestep

          rsub[j]       = frz1[i].r;
          uout[j]       = frz1[i].u;
          bulksub[j]    = frz1[i].bulk;
          shearsub[j]   = frz1[i].shear;
          shear33sub[j] = frz1[i].shear33;

          gradPsub      = frz1[i].gradP;
          inside        = frz1[i].inside;
          sigsub        = frz1[i].sigma;
          thetasub      = frz1[i].theta;
          Tfluc[j]      = frz1[i].T;             // replace with e
          muBfluc[j]    = frz1[i].muB;
          muSfluc[j]    = frz1[i].muS;
          muQfluc[j]    = frz1[i].muQ;
          sFO[j]        = frz1[i].s;
          eosname[j]    = frz1[i].eos_name;
        }
        else if ( swit == 2 ) // if particle was closer to freeze-out two timesteps ago
        {
          if ( p.btrack != -1 )                   // if particle had neighbors
            tlist[j]    = taupp;                  // two timesteps ago
          else
            tlist[j]    = taupp - systemPtr->dt;  // otherwise, go back three timesteps

          rsub[j]       = frz2[i].r;
          uout[j]       = frz2[i].u;
          bulksub[j]    = frz2[i].bulk;
          shearsub[j]   = frz2[i].shear;
          shear33sub[j] = frz2[i].shear33;

          gradPsub      = frz2[i].gradP;
          inside        = frz2[i].inside;
          sigsub        = frz2[i].sigma;
          thetasub      = frz2[i].theta;
          Tfluc[j]      = frz2[i].T;           // replace with e
          muBfluc[j]    = frz2[i].muB;
          muSfluc[j]    = frz2[i].muS;
          muQfluc[j]    = frz2[i].muQ;
          sFO[j]        = frz2[i].s;             // replace with e
          eosname[j]    = frz2[i].eos_name;
        }
        else  // this should never happen
        {
          cout << __PRETTY_FUNCTION__ << ": Not at freeze-out temperature" << endl;
        }

        // COMPUTE NORMALS AFTER THIS POINT
//        sFO[j]       = eosPtr->s_terms_T( Tfluc[j] );  // replace with e, BSQ
//cout << "Comparison: " << eosPtr->s_terms_T( Tfluc[j], muBfluc[j], muQfluc[j], muSfluc[j], eosname[j] )
//      << "   " << eosPtr->s_terms_T( Tfluc[j], eosname[j] ) << "   " << sFO[j] << endl;
        gsub[j]      = sqrt( Norm2(uout[j]) + 1 );


        sigsub      /= gsub[j]*tlist[j];
        swsub[j]     = p.norm_spec.s/sigsub;

        divT[j]      = (1.0/sFO[j])*gradPsub;
        divTtemp[j]  = -(1.0/(gsub[j]*sFO[j]))
                          *( cs2 * (wfz+bulksub[j]) * thetasub
                            - cs2*inside+inner(uout[j], gradPsub) );
    //THIS NEEDS TO BE RESET


        double insub = divTtemp[j]*divTtemp[j] - Norm2(divT[j]);
        double norm  = -sqrt(abs(insub));
        divTtemp[j] /= norm;
        divT[j]      = (1.0/norm)*divT[j];


        if ( divTtemp[j] == 1 )
        {
          cout << "track sph=" << p.btrack << " " << i << endl;
          cout << divTtemp[j] << " " << divT[j] << " " << norm << endl;
          cout << gradPsub << " " << thetasub << endl;
          cout << tlist[j] << " " << p.r << endl;
          cout << frz1[i].gradP << " " << frz2[i].gradP << endl;
          cout << frz1[i].T*197.3<< " " << frz2[i].T*197.3 << endl;
          getchar();
        }

        avgetasig += sFO[j]/sigsub;

        if(isnan(divTtemp[j]))
        {
          cout << "divtemp" << endl;
          cout << divTtemp[j] << " " << divT[j] << " " << norm << endl;
          cout << gradPsub << " " << thetasub << endl;
          cout << bulksub[j] << endl;
          cout << gsub[j] << endl;
          cout << tlist[j] << " " << p.r << endl;
          cout << frz1[i].T*0.1973<< " " << frz2[i].T*0.1973<< endl;
        }

        sFO[j]   *= pow(Tfluc[j]*0.1973, 3);
        Tfluc[j] *= 0.1973;

      }

      cf = curfrz;
    }


};

typedef std::shared_ptr<FreezeOut> pFreezeOut;  // smart pointer to freeze out object

#endif