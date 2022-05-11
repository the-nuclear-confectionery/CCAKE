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
  r(0)    = fields[0];
  r(1)    = fields[1];
  e_sub   = fields[2];
  rhoB_an = fields[3];
  rhoS_an = fields[4];
  rhoQ_an = fields[5];
  u(0)    = fields[6];
  u(1)    = fields[7];
  if ( fields.size() > 8 ) // passing in shear tensor initialization as well
  {
    double pi11 = fields[8];
    double pi22 = fields[9];
    double pi12 = fields[10];
    double pi33 = fields[11];
    shv(0,0)    = 0.0;
    shv(0,1)    = 0.0;
    shv(0,2)    = 0.0;
    shv(1,0)    = 0.0;
    shv(1,1)    = pi11;
    shv(1,2)    = pi12;
    shv(2,0)    = 0.0;
    shv(2,1)    = pi12;
    shv(2,2)    = pi22;
    shv33       = pi33;
  }
  s_an    = 0.0;
}

Particle::Particle( const Particle& p )
{
  Imat.identity();
  r       = p.r;
  e_sub   = p.e_sub;
  rhoB_an = p.rhoB_an;
  rhoS_an = p.rhoS_an;
  rhoQ_an = p.rhoQ_an;
  u       = p.u;
  shv     = p.shv;
  shv33   = p.shv33;
  s_an    = p.s_an;
}




////////////////////////////////////////////////////////////////////////////////
void Particle::set_SettingsPtr(Settings * settingsPtr_in)
{
  settingsPtr = settingsPtr_in;
}



////////////////////////////////////////////////////////////////////////////////
double Particle::gamcalc() { return sqrt( Norm2(u) + 1.0 ); }


////////////////////////////////////////////////////////////////////////////////
void Particle::frzcheck( double tin, int &count, int N )
{

// old version of freeze out
//  if ( Freeze == 0 )
//  {
//    if ( T() <= freezeoutT )
//    {
//      Freeze = 1;
//      frz2.t = tin;
//    }
//  }
//  else if ( Freeze == 1 )
//  {
//    if ( btrack == -1 )
//    {
//      count += 1;
//      Freeze = 3;
//      frz1.t = tin;
//    }
//    else if ( T()>frz1.T )
//    {
//      Freeze = 1;
//      frz2.t = tin;
//    }
//    else if( T() <= freezeoutT )
//    {
//      count += 1;
//      Freeze = 3;
//      frz1.t = tin;
//    }
//    else
//    {
//      Freeze=0;
//    }
//  }




// new version of freeze out

  if ( Freeze == 0 )
  {
    if ( e() <= efcheck )
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
    else if ( e()>frz1.e )
    {
      Freeze = 1;
      frz2.t = tin;
    }
    else if( e() <= efcheck )
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
void Particle::reset_pi_tensor(double tin2)
{
  gamma    = gamcalc();
  shv(2,1) = shv(1,2);
  shv(0,1) = 1./gamma*inner(u,colp1(1,shv));
  shv(0,2) = 1./gamma*inner(u,colp1(2,shv));
  shv(1,0) = shv(0,1);
  shv(2,0) = shv(0,2);

  mini(pimin,shv);
  uu       = u*u;
  shv(0,0) = 1./gamma/gamma*con(uu,pimin);
  shv33    = (shv(0,0)-shv(1,1)-shv(2,2))/tin2;
}




////////////////////////////////////////////////////////////////////////////////
// NEW FUNCTIONS TO CLEAN UP EQUATIONS_OF_MOTION BELOW THIS LINE
////////////////////////////////////////////////////////////////////////////////

void Particle::set_vartheta(double t) { vartheta = dsigma_dt/sigma - 1.0/t; }


////////////////////////////////////////////////////////////////////////////////
double Particle::get_Theta_force() { return -vartheta/gamma; }


////////////////////////////////////////////////////////////////////////////////
double Particle::get_Pi_force(double t)
{
  bigPI = Bulk*sigma/(gamma*t);
  return vartheta*(zeta/tauRelax + bigPI) - bigPI/(gamma*gamma*stauRelax);
}


////////////////////////////////////////////////////////////////////////////////
double Particle::get_aleph_force(double t)
{
  double aleph_force = -gamma*t*shv33;
  for (int i = 1; i < 3; i++)
  for (int j = 1; j < 3; j++)
    aleph_force += ( (shv(i,0)*u(j) + shv(0,j)*u(i))/gamma
                      - shv(i,j) - shv(0,0)*u(i)*u(j)/(gamma*gamma) )
                    * gradU(i,j);
  return aleph_force;
}


////////////////////////////////////////////////////////////////////////////////
Vector<double,2> Particle::get_Theta_mass()
{
  return v; // just the non-relativistic velocity
}


////////////////////////////////////////////////////////////////////////////////
Vector<double,2> Particle::get_Pi_mass()
{
  return (-(zeta/tauRelax + bigPI)/gamma)*v;
}


////////////////////////////////////////////////////////////////////////////////
Vector<double,2> Particle::get_aleph_mass()
{
  return shv(0,0)*v - colp1(0, shv);
}



////////////////////////////////////////////////////////////////////////////////
void Particle::set_hydro_info(double t)
{
  hydro.ID              = ID;
  hydro.t               = t;

  hydro.Agam            = Agam;
  hydro.Agam2           = Agam2;
  hydro.shv33           = shv33;

  hydro.gamma           = gamma;
  hydro.Bulk            = Bulk;
  hydro.bigPI           = bigPI;
  hydro.C               = C;
  hydro.tauRelax        = tauRelax;
  hydro.stauRelax       = stauRelax;
  hydro.zeta            = zeta;
  hydro.setas           = setas; 
  hydro.Ctot            = Ctot;
  hydro.Btot            = Btot;

  hydro.sigma           = sigma;
  hydro.dsigma_dt       = dsigma_dt;

  hydro.dwds            = dwds();
  hydro.dwdB            = dwdB();
  hydro.dwdS            = dwdS();
  hydro.dwdQ            = dwdQ();
  hydro.rhoB            = rhoB();
  hydro.rhoS            = rhoS();
  hydro.rhoQ            = rhoQ();

  hydro.g2              = g2;
  hydro.g3              = g3;
  hydro.gt              = gt;
  hydro.eta_o_tau       = eta_o_tau;
  hydro.dwdsT1          = dwdsT1;
  hydro.sigl            = sigl;

  hydro.T               = T();
  hydro.w               = w();
  hydro.s               = s();

  // vector members
  hydro.v               = v;
  hydro.u               = u;

  hydro.gradP           = gradP;
  hydro.gradBulk        = gradBulk;
  hydro.divshear        = divshear;
  hydro.gradshear       = gradshear;

  // matrix members
  hydro.Imat            = Imat;
  hydro.gradV           = gradV;
  hydro.gradU           = gradU;
  hydro.uu              = uu;
  hydro.pimin           = pimin;
  hydro.piu             = piu;
  hydro.piutot          = piutot;
  hydro.shv             = shv;
}



////////////////////////////////////////////////////////////////////////////////
void Particle::update_from_hydro_info()
{
  Agam = hydro.Agam;
  Agam2 = hydro.Agam2;

  bigPI = hydro.bigPI;
  C = hydro.C;
  Ctot = hydro.Ctot;
  Btot = hydro.Btot;

  sigma = hydro.sigma;
  dsigma_dt = hydro.dsigma_dt;

  g2 = hydro.g2;
  g3 = hydro.g3;
  gt = hydro.gt;
  eta_o_tau = hydro.eta_o_tau;
  dwdsT1 = hydro.dwdsT1;
  sigl = hydro.sigl;


  // vector members
//  v = hydro.v;
//  u = hydro.u;

  // matrix members
  gradU = hydro.gradU;
  uu = hydro.uu;
  pimin = hydro.pimin;
  piu = hydro.piu;
  piutot = hydro.piutot;

  div_u        = hydro.div_u;
  bigtheta     = hydro.bigtheta;
  inside       = hydro.inside;

  dBulk_dt     = hydro.dBulk_dt;
  detasigma_dt = hydro.detasigma_dt;
  du_dt        = hydro.du_dt;
  dshv_dt      = hydro.dshv_dt;
}