#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "particle.h"

////////////////////////////////////////////////////////////////////////////////
// Default constructor
Particle::Particle()
{
  s_an = 0.0;
  Imat.identity();
}
////////////////////////////////////////////////////////////////////////////////

// Overloaded constructor with initial fields
Particle::Particle(vector<double> &fields)
{
  Imat.identity();
  r.x[0]  = fields[0];
  r.x[1]  = fields[1];
  e_sub   = fields[2];
  rhoB_an = fields[3];
  rhoS_an = fields[4];
  rhoQ_an = fields[5];
  u.x[0]  = fields[6];
  u.x[1]  = fields[7];
  if ( fields.size() > 8 ) // passing in shear tensor initialization as well
  {
    double pi11 = fields[8];
    double pi22 = fields[9];
    double pi12 = fields[10];
    double pi33 = fields[11];
    shv.x[0][0] = 0.0;
    shv.x[0][1] = 0.0;
    shv.x[0][2] = 0.0;
    shv.x[1][0] = 0.0;
    shv.x[1][1] = pi11;
    shv.x[1][2] = pi12;
    shv.x[2][0] = 0.0;
    shv.x[2][1] = pi12;
    shv.x[2][2] = pi22;
    shv33       = pi33;
  }
  s_an    = 0.0;
}

Particle::Particle( const Particle& p )
{
  Imat.identity();
  r.x[0]  = p.r.x[0];
  r.x[1]  = p.r.x[1];
  e_sub   = p.e_sub;
  rhoB_an = p.rhoB_an;
  rhoS_an = p.rhoS_an;
  rhoQ_an = p.rhoQ_an;
  u.x[0]  = p.u.x[0];
  u.x[1]  = p.u.x[1];
  shv     = p.shv;
  shv33   = p.shv33;
  s_an    = p.s_an;
  eosPtr  = p.eosPtr;
}



////////////////////////////////////////////////////////////////////////////////
void Particle::set_EquationOfStatePtr( EquationOfState * eosPtr_in )
{
  eosPtr = eosPtr_in;
}

////////////////////////////////////////////////////////////////////////////////
void Particle::set_SettingsPtr(Settings * settingsPtr_in)
{
  settingsPtr = settingsPtr_in;
}

////////////////////////////////////////////////////////////////////////////////
double Particle::locate_phase_diagram_point_eBSQ(// previously s_out
                 double e_In, double rhoB_In, double rhoS_In, double rhoQ_In )
{
  double current_sVal = s();
  if ( Freeze == 5 )
    return current_sVal;
  else
  {
    // default: use particle's current location as initial guess
    eosPtr->tbqs( thermo.T, thermo.muB, thermo.muQ, thermo.muS );
    double sVal = eosPtr->s_out( e_In, rhoB_In, rhoS_In, rhoQ_In );

    if ( sVal > 0.0 )
      thermo.set(*eosPtr);
    else
      Freeze = 5; // new label for (totally decoupled) particles which go outside grid

    return sVal;
  }
}


////////////////////////////////////////////////////////////////////////////////
double Particle::locate_phase_diagram_point_eBSQ(double e_In)// previously s_out
{
  double current_sVal = s();
  if ( Freeze == 5 )
    return current_sVal;
  else
  {
    // default: use particle's current location as initial guess
    eosPtr->tbqs( thermo.T, 0.0, 0.0, 0.0 );
    double sVal = eosPtr->s_out( e_In, 0.0, 0.0, 0.0 );

    if ( sVal > 0.0 )
      thermo.set(*eosPtr);
    else
      Freeze = 5; // new label for (totally decoupled) particles which go outside grid

    return sVal;
  }
}

////////////////////////////////////////////////////////////////////////////////
void Particle::locate_phase_diagram_point_sBSQ(// previously update_s
                 double s_In, double rhoB_In, double rhoS_In, double rhoQ_In )
{
  if ( Freeze != 5 )
  {
    if ( abs(r.x[0]) < 0.00001 && abs(r.x[1]) < 0.00001 )
    cout << "CHECK SBSQ(1): " << s_In << "   " << rhoB_In << "   " << rhoS_In
        << "   " << rhoQ_In << endl;
    // default: use particle's current location as initial guess
    eosPtr->tbqs( thermo.T, thermo.muB, thermo.muQ, thermo.muS );
    bool update_s_success = eosPtr->update_s( s_In, rhoB_In, rhoS_In, rhoQ_In );

    if ( update_s_success )
      thermo.set(*eosPtr);
    else
      Freeze = 5; // new label for (totally decoupled) particles which go outside grid

    if ( abs(r.x[0]) < 0.00001 && abs(r.x[1]) < 0.00001 )
    cout << "CHECK SBSQ(2): " << s_In << "   " << rhoB_In << "   " << rhoS_In
        << "   " << rhoQ_In << endl;
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////
void Particle::locate_phase_diagram_point_sBSQ(double s_In) // previously update_s
{
  if ( Freeze != 5 )
  {
    // default: use particle's current location as initial guess
    eosPtr->tbqs( thermo.T, 0.0, 0.0, 0.0 );
    bool update_s_success = eosPtr->update_s( s_In, 0.0, 0.0, 0.0 );

    if ( update_s_success )
      thermo.set(*eosPtr);
    else
      Freeze = 5; // new label for (totally decoupled) particles which go outside grid
  }
  return;
}



////////////////////////////////////////////////////////////////////////////////
double Particle::gamcalc()
{
    return sqrt( Norm2(u) + 1.0 );
}


////////////////////////////////////////////////////////////////////////////////
void Particle::frzcheck( double tin, int &count, int N )
{
  if ( Freeze == 5 )
  {
    count++;
    return;
  }
  else if ( Freeze == 0 )
  {
    if ( T() <= freezeoutT )
    {
      Freeze = 1;
      frz2.t = tin;
    }
  }
  else if ( Freeze == 1 )
  {
    if ( btrack == -1 )
    {
      count += 1;
      Freeze = 3;
      frz1.t = tin;
    }
    else if ( T()>frz1.T )
    {
      Freeze = 1;
      frz2.t = tin;
    }
    else if( T() <= freezeoutT )
    {
      count += 1;
      Freeze = 3;
      frz1.t = tin;
    }
    else
    {
      Freeze=0;
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
//  Computes gamma and velocity
void Particle::calcbsq(double tin)
{
  gamma           = gamcalc();
  v               = (1.0/gamma)*u;
  double s_lab    = eta/gamma/tin;
//  double rhoB_lab = rhoB_sub/gamma/tin;
//  double rhoS_lab = rhoS_sub/gamma/tin;
//  double rhoQ_lab = rhoQ_sub/gamma/tin;
  double rhoB_lab = rhoB_sub;
  double rhoS_lab = rhoS_sub;
  double rhoQ_lab = rhoQ_sub;
  qmom            = ( (e()+p())*gamma/sigma )*u;
/*//THIS IS THE OLD ROUGH VERSION, JUST USE THE EVOLUTION EQUATIONS INSTEAD
	double rhoB_in2 = B*sigma/(gamma*sigmaweight);
	double rhoS_in2 = S*sigma/(gamma*sigmaweight);
	double rhoQ_in2 = Q*sigma/(gamma*sigmaweight);
	rhoB_an = rhoB_in2;
	rhoS_an = rhoS_in2;
	rhoQ_an = rhoQ_in2;
//cout << "\t - finding EoS solution for sBSQ: " << tin << "   " << s_lab << "   "
//		<< rhoB_in2 << "   " << rhoS_in2 << "   " << rhoQ_in2 << endl;
	locate_phase_diagram_point_sBSQ( s_lab, rhoB_in2, rhoS_in2, rhoQ_in2 );
*/
  // using *_sub quantities (which have been smoothed) in order to be consistent
  // with s_in2 = eta/(gamma*t) and evaluation of rest of EoM's quantities (e.g.,
  // sigma) which are also smoothed
	locate_phase_diagram_point_sBSQ( s_lab, rhoB_lab, rhoS_lab, rhoQ_lab );
}



////////////////////////////////////////////////////////////////////////////////
void Particle::return_bsqsv_A()
{
    eta_o_tau = (settingsPtr->using_shear) ? setas/stauRelax : 0.0;

	// THIS NEEDS TO BE CHECKED/FIXED, SPECIFICALLY WHEN INCLUDING MORE BETA-DOT TERMS
    Agam  = w() - dwds()*(s()+ bigPI/T() )- zeta/tauRelax
            - dwdB() * rhoB() - dwdS() * rhoS() - dwdQ() * rhoQ();

    //Agam2 = ( Agam - eta_o_tau*(0.5-1.0/3.0) - dwdsT1*shv.x[0][0] ) / gamma;
    Agam2 = ( Agam - eta_o_tau/6.0 - dwdsT1*shv.x[0][0] ) / gamma;
    Ctot  = C + 0.5*eta_o_tau*(1/g2-1);

}



////////////////////////////////////////////////////////////////////////////////
double Particle::Bsub()
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
    double pig  = shv.x[0][0]/g2;

    for (int i=0; i<=1; i++)
    for (int j=0; j<=1; j++)
      bsub += gradU.x[i][j] * ( pimin.x[i][j] + pig*uu.x[j][i]
                                - ( piu.x[i][j] + piu.x[j][i] ) / gamma );
    return bsub;
  }
}

////////////////////////////////////////////////////////////////////////////////
Matrix<double,2,2> Particle::Msub(int i)
{
  piu                     = rowp1(0,shv)*u;
  Matrix<double,2,2> msub = Agam2*uu + Ctot*gamma*Imat - (1+4/3./g2)*piu
                            + dwdsT1*transpose(piu) + gamma*pimin;
  return msub;
}


////////////////////////////////////////////////////////////////////////////////
Matrix<double,2,2> Particle::dpidtsub()
{
  Matrix<double,2,2> vsub;

  for (int i=0; i<=1; i++)
  for (int j=0; j<=1; j++)
  for (int k=0; k<=1; k++)
    vsub.x[i][j] += ( u.x[i]*pimin.x[j][k] + u.x[j]*pimin.x[i][k] )*du_dt.x[k];

  return vsub;
}

////////////////////////////////////////////////////////////////////////////////
void Particle::bsqsvsigset( double tin, int i )
{
  // from svsigset
  g2           = gamma*gamma;
  g3           = gamma*g2;
  gt           = gamma*tin;
  double dwdsT = dwds()/T();
  dwdsT1       = 1 - dwds()/T();
  sigl         = dsigma_dt/sigma - 1/tin;
  gradU        = gamma*gradV+g3*(v*(v*gradV));
  bigPI        = Bulk*sigma/gt ;
  C            = w()+ bigPI;
  return_bsqsv_A();
  Btot         = ( Agam*gamma + eta_o_tau/3*gamma )*sigl
                + bigPI/tauRelax + dwdsT*( gt*shv33 + Bsub() );
  check        = sigl;
}

////////////////////////////////////////////////////////////////////////////////
void Particle::setvisc( int etaconst, double bvf, double svf, double zTc,
                        double sTc, double sig, int type )
{
  //cout << __FUNCTION__ << ": " << etaconst << "   " << bvf << "   " << svf
  //      << "   " << zTc << "   " << sTc << "   " << sig << "   " << type << endl;

  if (type==1) // bulk viscosity
  {
    double temp=T()*197.3;
    zeta = bvf/(sig*sqrt(2*PI))*exp(-pow(temp-zTc,2)/(2.*sig*sig));
    zeta *= s();
    if (zeta<0.001) zeta=0.001;
    tauRelax = 9*zeta/(e()-3*p());
    if (tauRelax < 0.1) tauRelax=0.1;
  }
  else if (type==2) // shear viscosity
  {
    setas = svf*0.08;
    stauRelax=5*setas/w();
  }
  else if (type==3) // shear+bulk viscosity
  {
    if ((etaconst==1)||(etaconst==3)||(etaconst==4))
    {   // const eta/s
      setas=2*s()*svf;  // svf defines eta/s const (the two is needed for the
                           // definition in the code, don't remove!
    }
    //    for TECHQM/Gubser set svf=0.08
  }
  else if (type==4) // BSQ+shear+bulk viscosity
  {
    if ((etaconst==1)||(etaconst==3)||(etaconst==4))
    { // const eta/s
      setas=2*s()*svf;  // svf defines eta/s const; the two is needed
      // for the definition in the code, don't remove!
      //    for TECHQM/Gubser set svf=0.08
    }
    else
    { // eta/s (T)
      if (etaconst==5)
      {
        double TC=173.9; // 173.9/197.3
        double temp=T()*197.3/TC;
        double TC0=149.4/TC; // 149.4/197.3
        if( temp<TC0)
        {
          setas = s()*(8.0191 - 16.4659* temp  +  8.60918* temp *temp  );
        }
        else if( temp>1)
        {
          //setas = s()*(0.48 - 0.36*temp );
          setas = s()*(0.007407515123054544 +  0.06314680923610914* temp
                          + 0.08624567564083624* temp *temp  );
        }
        else
        {
          //setas = s()*(-0.107143 + 0.227143*temp);
          setas = s()*(0.397807 + 0.0776319* temp - 0.321513* temp *temp  );
        }
      }
      if (etaconst==6)
      {
        double TC=155; // 173.9/197.3
        double temp=T()*197.3/TC;
        double z=pow(0.66*temp,2);
        double alpha=33./(12.*PI)*(z-1)/(z*log(z));

        setas = s()*(0.0416762/pow(alpha,1.6)+ 0.0388977/pow(temp,5.1) );
      }
      else
      {
        double TC=sTc/197.3;
        double temp=T()/TC;
        if( temp>TC )
        {
          setas = s()*(0.3153036437246963 + 0.051740890688259315* temp
                          - 0.24704453441295576* temp *temp  );
        }
        else
        {
          setas = s()*(0.0054395278010882795 + 0.08078575671572835*temp
                          + 0.033774715483183566* temp *temp );
        }
      }
    }
    stauRelax=5*setas/w();
    if (!settingsPtr->using_Gubser && stauRelax < 0.005) stauRelax = 0.005;

    /// defining bulk viscosity

    if (bvf==0)
    {
      zeta =0;
      tauRelax = 1;
    }
    else
    {
      double temp=T()*197.3;

      if ((etaconst==2)||(etaconst==3))
      {
        double t2=temp/zTc;
        double min1=t2-1;
        if (t2>1.05) zeta=0.9*exp(-min1/0.025)+0.25*exp(-min1/0.13)+0.001;
        else if (t2<0.995) zeta=0.9*exp(min1/0.0025)+0.22*exp(min1/0.022)+0.03;
        else zeta=-13.77*t2*t2+27.55*t2-13.45;

        // single-argument version of cs2out
        tauRelax =5.*zeta/(pow((1-cs2()),2)*(e()+p()));
      }
      else if (etaconst==4)
      {
        double t2=temp/zTc;
        zeta=0.01162/sqrt(pow((t2-1.104),2)+ 0.0569777  ) - 0.1081/(t2*t2+23.7169);
        tauRelax = 5*( -0.0506*sTc/(temp*temp)
                      + 10.453/( temp*sqrt(0.156658 + pow( (temp/sTc-1.131),2) ) ) );
      }
      else
      {
        double t2=temp/zTc;
        double min1=t2-1;
        if (t2>1.05) zeta=0.9*exp(-min1/0.025)+0.25*exp(-min1/0.13)+0.001;
        else if (t2<0.995) zeta=0.9*exp(min1/0.0025)+0.22*exp(min1/0.022)+0.03;
        else zeta=-13.77*t2*t2+27.55*t2-13.45;

        // single-argument version of cs2out
        tauRelax =5.*zeta/(pow((1-cs2()),2)*(e()+p()));
      }

      zeta*=s();
      if (zeta<0.001) zeta=0.001;
      if (tauRelax <0.2) tauRelax=0.2;
    }

  }
}



////////////////////////////////////////////////////////////////////////////////
void Particle::sets(double tin2, bool is_first_timestep)
{
    gamma = gamcalc();
    if (settingsPtr->using_Gubser_with_shear && is_first_timestep)
    {
      shv.x[0][1] = 1./gamma*inner(u,colp1(1,shv));
      shv.x[0][2] = 1./gamma*inner(u,colp1(2,shv));
      shv.x[1][0] = shv.x[0][1];
      shv.x[2][0] = shv.x[0][2];
      
      setvar();
      //cout << "Sanity check: " << 1./gamma/gamma*con(uu,pimin) << "   "
      //      << shv.x[1][1] + shv.x[2][2] + tin2*shv33 << endl;

      //shv.x[0][0] = shv.x[1][1] + shv.x[2][2] + tin2*shv33;
      
      // go back to Jaki's default
      shv.x[0][0]=1./gamma/gamma*con(uu,pimin);
      shv33=(shv.x[0][0]-shv.x[1][1]-shv.x[2][2])/tin2;
    }
    else
    {
      shv.x[2][1]=shv.x[1][2];
      shv.x[0][1]=1./gamma*inner(u,colp1(1,shv));
      shv.x[0][2]=1./gamma*inner(u,colp1(2,shv));
      shv.x[1][0]=shv.x[0][1];
      shv.x[2][0]=shv.x[0][2];

      setvar();
      shv.x[0][0]=1./gamma/gamma*con(uu,pimin);
      shv33=(shv.x[0][0]-shv.x[1][1]-shv.x[2][2])/tin2;
    }
}


////////////////////////////////////////////////////////////////////////////////
void Particle::setvar()
{
    mini(pimin,shv);
    uu=u*u;
}

