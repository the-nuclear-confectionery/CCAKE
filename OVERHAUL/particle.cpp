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
    // default: use particle's current location as initial guess (pass in corresponding EoS as well!)
    eosPtr->tbqs( thermo.T, thermo.muB, thermo.muQ, thermo.muS, thermo.eos_name );
    bool solution_found = false;
    double sVal = eosPtr->s_out( e_In, rhoB_In, rhoS_In, rhoQ_In, solution_found );

    // save results if either default or conformal EoS returned a result
    // (assumes latter always succeeds)
    if ( solution_found )
      thermo.set(*eosPtr);
    else
      Freeze = 5; // new label for (totally decoupled) particles which go outside grid

    return sVal;
  }
}

////////////////////////////////////////////////////////////////////////////////
double Particle::locate_phase_diagram_point_eBSQ(double e_In)// previously s_out
                 { return locate_phase_diagram_point_eBSQ( e_In, 0.0, 0.0, 0.0 ); }






////////////////////////////////////////////////////////////////////////////////
void Particle::locate_phase_diagram_point_sBSQ(// previously update_s
                 double s_In, double rhoB_In, double rhoS_In, double rhoQ_In )
{
  if ( Freeze != 5 )
  {
    // default: use particle's current location as initial guess
    eosPtr->tbqs( thermo.T, thermo.muB, thermo.muQ, thermo.muS, thermo.eos_name );
    bool update_s_success = eosPtr->update_s( s_In, rhoB_In, rhoS_In, rhoQ_In );

    if ( update_s_success )
      thermo.set(*eosPtr);
    else
      Freeze = 5; // new label for (totally decoupled) particles which go outside grid
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////
void Particle::locate_phase_diagram_point_sBSQ(double s_In) // previously update_s
               { locate_phase_diagram_point_sBSQ( s_In, 0.0, 0.0, 0.0 ); }



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
//  Computes gamma and velocity
void Particle::calcbsq(double tin)
{
  gamma           = gamcalc();
  v               = (1.0/gamma)*u;
  double s_lab    = eta/gamma/tin;
  double rhoB_lab = rhoB_sub/gamma/tin;
  double rhoS_lab = rhoS_sub/gamma/tin;
  double rhoQ_lab = rhoQ_sub/gamma/tin;
  qmom            = ( (e()+p())*gamma/sigma )*u;
  // using *_sub quantities (which have been smoothed) in order to be consistent
  // with s_in2 = eta/(gamma*t) and evaluation of rest of EoM's quantities (e.g.,
  // sigma) which are also smoothed
//  std::cout << "Doing this particle: " << r << "   " << eta << std::endl;
	locate_phase_diagram_point_sBSQ( s_lab, rhoB_lab, rhoS_lab, rhoQ_lab );
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
    double pig  = shv(0,0)/g2;

    for (int i=0; i<=1; i++)
    for (int j=0; j<=1; j++)
      bsub += gradU(i,j) * ( pimin(i,j) + pig*uu(j,i) - ( piu(i,j) + piu(j,i) ) / gamma );
    return bsub;
  }
}

////////////////////////////////////////////////////////////////////////////////
Matrix<double,2,2> Particle::Msub()
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
    vsub(i,j) += ( u(i)*pimin(j,k) + u(j)*pimin(i,k) )*du_dt(k);

  return vsub;
}

////////////////////////////////////////////////////////////////////////////////
void Particle::update_fluid_quantities( double tin )
{
  // from svsigset
  g2           = gamma*gamma;
  g3           = gamma*g2;
  gt           = gamma*tin;
  double dwdsT = dwds()/T();
  dwdsT1       = 1 - dwds()/T();
  sigl         = dsigma_dt/sigma - 1/tin;
  gradU        = gamma*gradV+g3*(v*(v*gradV));
  bigPI        = Bulk*sigma/gt;
  C            = w()+ bigPI;

  eta_o_tau    = (settingsPtr->using_shear) ? setas/stauRelax : 0.0;

  // THIS NEEDS TO BE CHECKED/FIXED,
  // SPECIFICALLY WHEN INCLUDING MORE BETA-DOT TERMS
  Agam         = w() - dwds()*( s()+ bigPI/T() )- zeta/tauRelax
                 - dwdB() * rhoB() - dwdS() * rhoS() - dwdQ() * rhoQ();

  Agam2        = ( Agam - eta_o_tau/3.0 - dwdsT1*shv(0,0) ) / gamma;
  Ctot         = C + eta_o_tau*(1.0/g2-1.0);


  Btot         = ( Agam*gamma + 2.0*eta_o_tau/3.0*gamma )*sigl
                + bigPI/tauRelax + dwdsT*( gt*shv33 + Bsub() );
  check        = sigl;
//std::cout << "CHECK update_fluid_quantities: " << ID << "   " << tin << "   " << g2
//          << "   " << g2 << "   " << g3 << "   " << gt << "   "
//          << dwds() << "   " << T() << "   " << sigl << "   "
//          << w() << "   " << bigPI << "   " << gradU << "   "
//          << Agam << "   " << eta_o_tau << "   " << tauRelax << "   "
//          << shv33 << "   " << Bsub() << "   " << Btot << "   "
//          << Agam2 << "   " << Ctot << "\n";
}

////////////////////////////////////////////////////////////////////////////////
void Particle::setvisc( int etaconst, double bvf, double svf, double zTc,
                        double sTc, double sig, int type )
{
  setas     = s()*svf;

  stauRelax = 5*setas/w();
  if (!settingsPtr->using_Gubser && stauRelax < 0.005) stauRelax = 0.005;

  /*if (abs(bvf) > 1e-6)
  {
    cerr << "You need to replace setvisc!!!" << endl;
    exit(1);
  }
  zeta = 0;
  tauRelax = 1;*/
  zeta = s()*bvf;
  tauRelax = 5.*zeta/(pow((1.0-cs2()),2.0)*w());
  if (tauRelax < 0.1) tauRelax=0.1;
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




void Particle::evaluate_time_derivatives( double t )
{

  double gamt = 0.0, pre = 0.0, p1 = 0.0;
  if ( settingsPtr->using_shear )
  {
    gamt = 1.0/gamma/stauRelax;
    pre  = eta_o_tau/gamma;
    p1   = gamt - 4.0/3.0/sigma*dsigma_dt + 1.0/t/3.0;
  }

  Vector<double,2> minshv   = rowp1(0, shv);
  Matrix <double,2,2> partU = gradU + transpose( gradU );

  //===============
  // print status
  if ( VERBOSE > 2 && print_this_particle )
  {
    std::cout << "CHECK misc1: " << ID << "   " << t << "   "
              << gamt << "   " << sigma	<< "   " << dsigma_dt << "\n"
              << "CHECK minshv: " << ID << "   " << t << "   " << minshv << "\n"
              << "CHECK partU: " << ID << "   " << t << "   " << partU << "\n";
  }


  // set the Mass and the Force
  Matrix <double,2,2> M = Msub();
  Vector<double,2> F    = Btot*u + gradshear - ( gradP + gradBulk + divshear );


  //===============
  // print status
  if ( VERBOSE > 2 && print_this_particle )
    std::cout << "CHECK M: " << ID << "   " << t << "   " << M << "\n"
              << "CHECK F: " << ID << "   " << t << "   " << F << "   "
              << Btot << "   " << u << "   "
              << gradshear << "   " << gradP << "   "
              << gradBulk << "   " << divshear << "\n";

  //===============
  // shear contribution
  if ( settingsPtr->using_shear )
    F += pre*v*partU + p1*minshv;

  //===============
  // print status
  if ( VERBOSE > 2 && print_this_particle )
    std::cout << "CHECK F(again): " << ID << "   " << t << "   "
              << F << "   " << pre << "   " << v << "   " << partU << "   "
              << p1 << "   " << minshv << "\n";


  double det = deter(M);

  Matrix <double,2,2> MI;
  MI(0,0) =  M(1,1)/det;
  MI(0,1) = -M(0,1)/det;
  MI(1,0) = -M(1,0)/det;
  MI(1,1) =  M(0,0)/det;


  //===============
  // print status
  if ( VERBOSE > 2 && print_this_particle )
    std::cout << "CHECK det: " << ID << "   " << t << "   "
              << M << "   " << det << "\n"
              << "CHECK MI: " << ID << "   " << t
              << "   " << MI << "\n";


  //===============
  // compute acceleration
  du_dt(0) = F(0) * MI(0,0) + F(1) * MI(0,1);
  du_dt(1) = F(0) * MI(1,0) + F(1) * MI(1,1);

  //===============
  // define auxiliary variables
  double vduk               = inner( v, du_dt );
  Matrix <double,2,2> ulpi  = u*colp1(0, shv);
  Matrix <double,2,2> Ipi   = -2.0*eta_o_tau/3.0 * ( Imat + uu ) + 4./3.*pimin;

  //===============
  // "coordinate" divergence
  div_u                   = (1./ gamma)*inner( u, du_dt)
                                - ( gamma/ sigma ) * dsigma_dt;
  //===============
  // "covariant" divergence
  bigtheta                = div_u*t+gamma;
      /* the above lines could automatically be set in particle after 
      calculating the matrix elements above */

  //===============
  // print status
  if ( VERBOSE > 2 && print_this_particle )
    std::cout << "CHECK div_u: " << ID
              << "   " << t
              << "   " << div_u
              << "   " << gamma
              << "   " << u
              << "   " << du_dt
              << "   " << inner( u, du_dt)
              << "   " << sigma 
              << "   " << dsigma_dt << "\n"
              << "CHECK bigtheta: " << ID
              << "   " << t
              << "   " << bigtheta
              << "   " << gamma << "\n";

  //===============
  // this term occurs in Eqs. (250) and (251) of Jaki's long notes
  // translation: pi^{ij} + pi^{00} v^i v^j - pi^{i0} v^j - pi^{0j} v^i
  Matrix <double,2,2> sub   = pimin + (shv(0,0)/g2)*uu -1./gamma*piutot;

  // minshv = pi^{0i}                   (i   = 1,2)
  // pimin  = pi^{ij}                   (i,j = 1,2)
  // uu     = u^i u^j                   (i,j = 1,2)
  // piu    = pi^{0i} u^j               (i,j = 1,2)
  // piutot = pi^{0i} u^j + pi^{0j} u^i (i,j = 1,2)
  // gradU  = du_i/dx^j                 (i,j = 1,2)

  //===============
  if ( settingsPtr->using_shear )
    inside                  = t*( inner( -minshv+shv(0,0)*v, du_dt )
                                  - con2(sub, gradU) - gamma*t*shv33 );


  //===============
  // print status
  if ( VERBOSE > 2 && print_this_particle )
    std::cout << "CHECK inside: " << ID << "   "
              << t << "   "
              << inside << "   "
              << minshv << ";   "
              << shv(0,0)*v << ";   "
              << du_dt << ";   "
              << sub << "   "
              << gradU << ";   "
              << gamma*t*shv33 << "\n";



  detasigma_dt            = 1./sigma/T()*( -bigPI*bigtheta + inside );


  //===============
  // print status
  if ( VERBOSE > 2 && print_this_particle )
    std::cout << "CHECK detasigma_dt: " << ID << "   "
              << t << "   "
              << detasigma_dt << "   "
              << sigma << "   "
              << T()*hbarc_MeVfm << "   "
              << bigPI << "   "
              << bigtheta << "   "
              << inside << "\n";


  // N.B. - ADD EXTRA TERMS FOR BULK EQUATION
  dBulk_dt = ( -zeta/sigma*bigtheta - Bulk/gamma )/tauRelax;

  Matrix <double,2,2> ududt = u*du_dt;

  // N.B. - ADD READABLE TERM NAMES
  if ( settingsPtr->using_shear )
    dshv_dt                 = - gamt*( pimin + setas*partU )
                             - eta_o_tau*( ududt + transpose(ududt) )
                             + dpidtsub() + sigl*Ipi
                             - vduk*( ulpi + transpose(ulpi) + (1/gamma)*Ipi );
}