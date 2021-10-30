#include "system_state.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <sstream>

//using namespace std;
using std::cout;
using std::endl;
using std::string;

#include "constants.h"
#include "vector.h"
#include "particle.h"
#include "runge_kutta.h"
#include "eos.h"
#include "Stopwatch.h"
#include "system_state.h"
#include "linklist.h"

using namespace constants;


////////////////////////////////////////////////////////////////////////////////
void SystemState::set_EquationOfStatePtr( EquationOfState * eosPtr_in )
{
  eosPtr = eosPtr_in;
}


void SystemState::set_SettingsPtr(Settings * settingsPtr_in)
{
  settingsPtr = settingsPtr_in;
}


////////////////////////////////////////////////////////////////////////////////
void SystemState::initialize()  // formerly called "manualenter"
{
  int start, end;
  int df;

  cout << __PRETTY_FUNCTION__ << ": " << endl;
  _n = particles.size();

  cout << "_n = " << _n << endl;

  t = settingsPtr->t0;

  cout << "t = " << t << endl;

  _h = settingsPtr->_h;

  settingsPtr->efcheck = eosPtr->efreeze(settingsPtr->Freeze_Out_Temperature);
  settingsPtr->sfcheck = eosPtr->sfreeze(settingsPtr->Freeze_Out_Temperature);

		std::cout << "efcheck = " << settingsPtr->efcheck*hbarc_GeVfm << " GeV/fm^3\n";
		std::cout << "sfcheck = " << settingsPtr->sfcheck << " 1/fm^3\n";


  for (auto & p : particles) p.set_EquationOfStatePtr( eosPtr );

  linklist.efcheck = efcheck;
  linklist.sfcheck = sfcheck;
  linklist.fcount  = 0;
  linklist.average = 0;
  //       Start reading ICs          //

  //int numpart, _Ntable3;

  //  cout << "setting up SPH" << endl;
  return;
}



void SystemState::initialize_linklist()
{

  string ictype = "iccing";
  cout << "Initial conditions type: " << ictype << endl;

  if ( ictype == "iccing" )
  {
    settingsPtr->gtyp=6;

    int count           = 1;
    vector<string>        filelist( count );

    int j               = 0;
    filelist[j]         = "./ic0.dat"; // only doing single event
    linklist.filenames  = filelist;
    linklist.fcount     = count;
    linklist.fnum       = linklist.start;
    
    cout << "Check 0: " << particles[0].r.x[0] << "   " << particles[0].r.x[1] << endl;

    int currently_frozen_out = number_part;
    linklist.initialize( settingsPtr->t0, particles.size(),
                         settingsPtr->_h, &particles, dt, currently_frozen_out );

    //cout << "number of sph particles=" << _Ntable3 << endl;
    linklist.gtyp=settingsPtr->gtyp;

  }


  // formerly bsqsv_set in this loop
  for (auto & p : particles)
  {
    double gg = p.gamcalc();
    p.g2      = gg*gg;
    p.shv33   = 0.0;
  }

  return;
}






///////////////////////////////////////////////////////////////////////////////
//Dekra: BSQsimulation was moved from this location to BSQhydro and renamed to run
///////////////////////////////////////////////////////////////////////////////
//Start routines for checking conservation of energy density, entropy density and BSQ 
// charge densities at each time step of the simulation
///////////////////////////////////////
void SystemState::check_BSQ_energy_conservation()
{
  E=0.0;
  for ( auto & p : particles )
    E += ( p.C*p.g2 - p.p() - p.bigPI + p.shv.x[0][0] )
          *p.sigmaweight*t/p.sigma;

  if (linklist.first == 1)
  {
    linklist.first = 0;
    E0    = E;
  }

  return;
}
////////////////////////////////////////
void SystemState::check_BSQ_charge_conservation()
{
  Btotal = 0.0;
  Stotal = 0.0;
  Qtotal = 0.0;

  for ( auto & p : particles )
  {
    //Btotal += p.B;
    //Stotal += p.S;
    //Qtotal += p.Q;
    Btotal += p.rhoB_sub*p.rho_weight;
    Stotal += p.rhoS_sub*p.rho_weight;
    Qtotal += p.rhoQ_sub*p.rho_weight;
  }

  if (linklist.first==1)
  {
    Btotal0 = Btotal;
    Stotal0 = Stotal;
    Qtotal0 = Qtotal;
  }

	return;
}
///////////////////////////////////////
void SystemState::bsqsvconservation()
{
    bsqsvconservation_E();
    Etot  = E + Ez;
    Eloss = (E0-Etot)/E0*100;
    rk2   = 0;
}
///////////////////////////////////////
void SystemState::conservation_entropy()
{
  S=0.0;

  for (int i=0; i<_n; i++)
  {
    S += particles[i].eta_sigma*particles[i].sigmaweight;
    if (i==0)
    std::cout << "\t\t --> " << i << "   " << particles[i].eta_sigma << "   "
              << particles[i].sigmaweight << "   " << S << endl;
  }

  if (linklist.first==1)
    S0=S;
}
///////////////////////////////////////
void SystemState::conservation_BSQ()
{
    Btotal = 0.0;
    Stotal = 0.0;
    Qtotal = 0.0;

    for (int i=0; i<_n; i++)
	{
        //Btotal += particles[i].B;
        //Stotal += particles[i].S;
        //Qtotal += particles[i].Q;
        Btotal += particles[i].rhoB_sub*particles[i].rho_weight;
        Stotal += particles[i].rhoS_sub*particles[i].rho_weight;
        Qtotal += particles[i].rhoQ_sub*particles[i].rho_weight;
    }

    if (linklist.first==1)
    {
        Btotal0 = Btotal;
        Stotal0 = Stotal;
        Qtotal0 = Qtotal;
    }
	return;
}


///////////////////////////////////////
void SystemState::bsqsvconservation_E()
{

    E=0.;
    for (int i=0; i<_n; i++)
    {
      auto & p = particles[i];

        E += ( p.C*p.g2 - p.p() - p.bigPI + p.shv.x[0][0] )
              / p.sigma*p.sigmaweight*t;
        if (i==0)
          std::cout << "E: " << i << "   " << t
              << "   " << p.T()
              << "   " << p.e()
              << "   " << p.C
              << "   " << p.g2
              << "   " << p.p()
              << "   " << p.bigPI
              << "   " << p.shv.x[0][0]
              << "   " << p.sigma
              << "   " << p.sigmaweight << endl;    }

    if (linklist.first==1)
    {
      linklist.first=0;
      E0=E;
    }
}


///////////////////////////////////////
void SystemState::bsqsvconservation_Ez()
{
  dEz=0.;

  double t2=t*t;
  for (int i=0; i<_n; i++)
  {
    auto & p = particles[i];

    dEz += ( p.p() + p.bigPI + p.shv33*t2 ) / p.sigma*p.sigmaweight;

    if (false)
      std::cout << "dEz: " << i << "   " << t
        << "   " << p.p()
        << "   " << p.bigPI
        << "   " << p.shv33*t2
        << "   " << p.sigma
        << "   " << p.sigmaweight << endl;
  }
}






///////////////////////////////////////////////////////////////////////////////
void SystemState::set_current_timestep_quantities()
{
  N = _n;

  etasigma0.resize(N);
  Bulk0.resize(N);

  u0.resize(N);
  r0.resize(N);

  shv0.resize(N);

  cout << __PRETTY_FUNCTION__ << ": N = " << N << endl;

  for (int i=0; i<N; ++i)
  {
    const auto & p = particles[i];
    u0[i]        = p.u;
    r0[i]        = p.r;
    etasigma0[i] = p.eta_sigma;
    Bulk0[i]     = p.Bulk;
    mini( shv0[i], p.shv );
  }
}
///////////////////////////////////////////////////////////////////////////////
void SystemState::get_derivative_halfstep(double dx)
{
  N = _n;

  for (int i=0; i<N; ++i)
  {
    auto & p = particles[i];
    p.u            = u0[i]        + 0.5*dx*p.du_dt;
    p.r            = r0[i]        + 0.5*dx*p.v;
    p.eta_sigma    = etasigma0[i] + 0.5*dx*p.detasigma_dt;
    p.Bulk         = Bulk0[i]     + 0.5*dx*p.dBulk_dt;
    tmini( p.shv,    shv0[i]      + 0.5*dx*p.dshv_dt );
  }
}
///////////////////////////////////////////////////////////////////////////////
void SystemState::get_derivative_fullstep(double dx)
{
  N = _n;

  for (int i=0; i<N; ++i)
  {
    auto & p = particles[i];
    p.u            = u0[i]        + dx*p.du_dt;
    p.r            = r0[i]        + dx*p.v;
    p.eta_sigma    = etasigma0[i] + dx*p.detasigma_dt;
    p.Bulk         = Bulk0[i]     + dx*p.dBulk_dt;
    tmini( p.shv,    shv0[i]      + dx*p.dshv_dt );
  }
}
//////////////////////////////////////////////////////////////////////////////



void SystemState::bsqsvfreezeout(int curfrz)
{
    if (frzc==0)
    {
        taupp=t;
        frzc=1;
        for (int i=0; i<_n; i++) {


            _p[i].frz2.r=_p[i].r;
            _p[i].frz2.u=_p[i].u;
            _p[i].frz2.sigma=_p[i].sigma;
            _p[i].frz2.T=_p[i].EOST();
            _p[i].frz2.bulk=_p[i].bigPI ;
            _p[i].frz2.theta=_p[i].div_u+_p[i].gamma/t;
            _p[i].frz2.gradP=_p[i].gradP;
            _p[i].frz2.shear=_p[i].shv;
            _p[i].frz2.shear33=_p[i].shv33;
            _p[i].frz2.inside=_p[i].inside;
        }

    }
    else if (frzc==1)
    {
        taup=t;
        frzc=2;
        for (int i=0; i<_n; i++) {

            _p[i].frz1.r=_p[i].r;
            _p[i].frz1.u=_p[i].u;
            _p[i].frz1.sigma=_p[i].sigma;
            _p[i].frz1.T=_p[i].EOST();
            _p[i].frz1.bulk=_p[i].bigPI ;
            _p[i].frz1.theta=_p[i].div_u+_p[i].gamma/t;
            _p[i].frz1.gradP=_p[i].gradP;
            _p[i].frz1.shear=_p[i].shv;
            _p[i].frz1.shear33=_p[i].shv33;
            _p[i].frz1.inside=_p[i].inside;
        }

        divTtemp=new double [curfrz];
        divT=new Vector<double,D> [curfrz];
        gsub=new double [curfrz];
        uout=new Vector<double,D> [curfrz];
        swsub=new double [curfrz];
        bulksub=new double [curfrz];
        shearsub=new Matrix<double,D+1,D+1> [curfrz];
        shear33sub=new double [curfrz];
        tlist=new double [curfrz];
        rsub=new Vector<double,D> [curfrz];

        if (curfrz>0)
            bsqsvinterpolate(curfrz);
        else
            cf=0;

    }
    else
    {

        for (int i=0; i<_n; i++) {
            if (_p[i].Freeze<4) {
                if ((_p[i].btrack<=3)&&(_p[i].btrack>0)) {
                    _p[i].fback4=_p[i].fback2;
                    _p[i].fback3=_p[i].fback;
                    _p[i].fback2=_p[i].frz2;
                    _p[i].fback=_p[i].frz1;
                }
                else if (_p[i].btrack==0) {
                    if (_p[i].fback.gradP.x[0]!=0) {
                        _p[i].frz2=_p[i].fback2;
                        _p[i].frz1=_p[i].fback;
                    }
                    else {
                        _p[i].frz2=_p[i].fback4;
                        _p[i].frz1=_p[i].fback3;
                        cout << "back second"  << endl;
                    }


                    curfrz++;
                    list.push_back(i);
                    _p[i].Freeze=4;
                    _p[i].btrack=-1;
                }
            }
        }

        tau=t;

        divTtemp=new double [curfrz];
        divT=new Vector<double,D> [curfrz];
        gsub=new double [curfrz];
        uout=new Vector<double,D> [curfrz];
        swsub=new double [curfrz];
        bulksub=new double [curfrz];
        shearsub=new Matrix<double,D+1,D+1> [curfrz];
        shear33sub=new double [curfrz];
        tlist=new double [curfrz];
        rsub=new Vector<double,D> [curfrz];



        if (curfrz>0)
            bsqsvinterpolate(curfrz);
        else
            cf=0;


        //sets up the variables for the next time step
        for (int i=0; i<_n; i++) {
            _p[i].frz2=_p[i].frz1;

            _p[i].frz1.r=_p[i].r;
            _p[i].frz1.u=_p[i].u;
            _p[i].frz1.sigma=_p[i].sigma;
            _p[i].frz1.T=_p[i].EOST();
            _p[i].frz1.bulk=_p[i].bigPI ;
            _p[i].frz1.theta=_p[i].div_u+_p[i].gamma/t;
            _p[i].frz1.gradP=_p[i].gradP;
            _p[i].frz1.shear=_p[i].shv;
            _p[i].frz1.shear33=_p[i].shv33;
            _p[i].frz1.inside=_p[i].inside;
        }
        taupp=taup;
        taup=tau;
    }
    cfon=0;
}




void SystemState::bsqsvinterpolate(int curfrz)
{

    sFO.resize(curfrz,0);
    Tfluc.resize(curfrz,0);
    for (int j=0; j<curfrz; j++)
    {
        int i=list[j];

        int swit=0;
        if (abs(_p[i].frz1.T-freezeoutT)<abs(_p[i].frz2.T-freezeoutT)) swit=1;
        else swit=2;

        double sigsub,thetasub,inside;
        Vector<double,D> gradPsub;
        if (swit==1) {
            if (_p[i].btrack!=-1) tlist[j]=taup;
            else tlist[j]=taup-dt;
            rsub[j]=_p[i].frz1.r;
            uout[j]=_p[i].frz1.u;
            bulksub[j]=_p[i].frz1.bulk;
            shearsub[j]=_p[i].frz1.shear;
            shear33sub[j]=_p[i].frz1.shear33;

            gradPsub=_p[i].frz1.gradP;
            inside=_p[i].frz1.inside;
            sigsub=_p[i].frz1.sigma;
            thetasub=_p[i].frz1.theta;
            Tfluc[j]=_p[i].frz1.T;
        }
        else if (swit==2) {
            if (_p[i].btrack!=-1) tlist[j]=taupp;
            else tlist[j]=taupp-dt;
            rsub[j]=_p[i].frz2.r;
            uout[j]=_p[i].frz2.u;
            bulksub[j]=_p[i].frz2.bulk;
            shearsub[j]=_p[i].frz2.shear;
            shear33sub[j]=_p[i].frz2.shear33;

            gradPsub=_p[i].frz2.gradP;
            inside=_p[i].frz2.inside;
            sigsub=_p[i].frz2.sigma;
            thetasub=_p[i].frz2.theta;
            Tfluc[j]=_p[i].frz2.T;
        }
        else {
            cout << "LinkList.h: Not at freeze-out temperature" << endl;

        }

        sFO[j]=_p[0].EOSs_terms_T(Tfluc[j]);

        gsub[j]=sqrt( Norm2(uout[j]) + 1 );


        sigsub/=gsub[j]*tlist[j];
        swsub[j]=_p[i].sigmaweight/sigsub;

        divT[j]=(1/sFO[j])*gradPsub;
        divTtemp[j]=-(1/(gsub[j]*sFO[j]))*(cs2*(wfz+bulksub[j])*thetasub-cs2*inside+inner(uout[j],gradPsub));


        double insub=divTtemp[j]*divTtemp[j]-Norm2(divT[j]);
        double norm=-sqrt(abs(insub));
        divTtemp[j]/=norm;
        divT[j]=(1/norm)*divT[j];



        if (divTtemp[j]==1) {
            cout << "track sph=" << _p[i].btrack << " " << i <<  endl;
            cout << divTtemp[j] << " " << divT[j] << " " << norm << endl;
            cout << gradPsub << " " << thetasub << endl;
            cout << tlist[j] << " " << _p[i].r << endl;
            cout << _p[i].frz1.gradP<< " " << _p[i].frz2.gradP <<  endl;
            cout << _p[i].frz1.T*197.3<< " " << _p[i].frz2.T*197.3 <<  endl;
            getchar();

        }

        avgetasig+=sFO[j]/sigsub;

        if(isnan(divTtemp[j]))
        {

            cout << "divtemp" << endl;
            cout << divTtemp[j] << " " << divT[j] << " " << norm << endl;
            cout << gradPsub << " " << thetasub << endl;
            cout << bulksub[j] <<endl;
            cout << gsub[j] << endl;
            cout << tlist[j] << " " << _p[i].r << endl;
            cout << _p[i].frz1.T*0.1973<< " " << _p[i].frz2.T*0.1973<<  endl;

        }

        sFO[j]*=pow(Tfluc[j]*0.1973,3);
        Tfluc[j]*=0.1973;

    }
    cf=curfrz;
}