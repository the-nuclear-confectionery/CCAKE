#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "vector.h"
#include "matrix.h"
#include "eos.h"
#include <stdio.h>
#include <math.h>

extern double freezeoutT;

// need to switch to dimension instead!
template <int D>
class Particle
{
private:


public:
////////////////////////////////////////////////////////////////////////////////
//                            Particle Variables                              //
////////////////////////////////////////////////////////////////////////////////


    int btrack;
    static eos EOS;	//use one copy of EOS for all particles
    double Agam, Agam2;
    double sigmaweight;        // specific volume per particle
    Vector<double,D> r;                   // position
    Vector<double,D> v;                   // velocity
    Vector<double,D> u;                   // relativistic velocity
    Vector<double,D> qmom;
    Vector<double,D> du_dt,gradsig;               // relativistic velocity derivative
    Vector<double,D> k1,k2,k3,k4,ksub;
    Vector<double,D> r1,r2,r3,r4,rsub;
    double ets1,ets2,ets3,ets4;
    double b1,b2,b3,b4;
    double bn1,bn2,bn3,bn4;
    Matrix<double,D,D> shv1,shv2,shv3,svh4;
    double shv33;
//    double vmag,vang;

    struct FRZ
    {
        Vector<double,D> r,u,gradP;
        double t,s,T,theta,bulk,sigma,shear33,inside;
        Matrix<double,D+1,D+1> shear;
    };

    FRZ frz1,frz2,fback,fback2,fback3,fback4;

//    double es1,es2,es3,es4;
    double div_u;              // four-divergence of relativistic velocity
    double gamma;               // Lorentz factor
    double gamcalc();
    //double gamcalcbsq();
    double s_sub,e_sub,s_an,s_rat,sigsub;
    double eta_sigma;          // Ratio entropy/especific volume
    double detasigma_dt;
    double Bulk;               // Bulk Viscosity weight
    double bigPI;        // total bulk viscosity
    double C;
    double tauRelax;           // Bulk Relaxation time
    double stauRelax;        // Shear Relxation time;
    double dBulk_dt;           // derivative Bulk Viscosity
    double zeta;               // bulk coefficient
    double setas;
    int count;
    int Freeze;
    double Ctot, Btot;
    Matrix<double,D+1,D+1> shv;
    Matrix<double,D,D> shv0,dshv_dt;
    Vector<double,D>  divshear,gradshear;

    double sv_eta,taupi;

////////////////////////////////////////////////////////////////////////////////
//                           Fluid Variables                                  //
////////////////////////////////////////////////////////////////////////////////
    //double sigmastar;
    double sigma;              // especific volume
    double dsigma_dt;          // derivative of especific volume
    //double div_J;              // divergence of the flux
    //Vector<double,D> flux;               // flux = sigma * velocity
    Vector<double,D> gradP;              // Gradient of Pressure
    Matrix<double,D,D> gradV,gradU;          // Gradient of velocity needed for shear

    Vector<double,D> gradBulk,gradrhoB,gradrhoS,gradrhoQ;           // Gradient of Bulk Viscosity
    Vector<double,D> gradsigma;          // Gradient of especific volume

    double dw_ds;              // derivative of the enthalpy on entropy
    double eta;                   // entropy density
    double rhoB;                   // Baryon density
    double rhoS;                   // strange density
    double rhoQ;                   // electric charge density
//    double P;                   // pressure
//    double epsilon;               // energy density
//    double s;                   // entropy density
//    double T;                   // temperature
    double eden;

    double bigtheta;

	double B, S, Q;						// baryon, strange, electric charge
    double drhoB_dt,drhoS_dt,drhoQ_dt;

    Particle();
    void start(string enter, eos & EOS_in);
    void calc( double tin );
    void calcbsq(double tin);
    void sigset(double tin);
    void vsigset(double tin);
    void svsigset(double tin, int i);
    void bsqsvsigset(double tin, int i);
    void returnv_A();
    void return_sv_A();
    void return_bsqsv_A();
    void returnA() ;
    void setvisc(int etaconst,double bvf, double svf, double zTc, double sTc, double sig, int type);
    void frzcheck(double tin,int &count, int N);
    double g2,g3,gt;
    double eta_o_tau,dwdsT1,sigl;
    Matrix <double,D,D> uu,pimin,piu,piutot;
    //void ucalc();
    double Bsub()  ;
    Matrix<double, D,D> Msub(int i),Imat  ;
    Matrix<double, D,D> dpidtsub();
    void sets(double tin2) ;
    double bcheck;
    double check;
    void setvar();
    double saves;
    double inside;

	double EOST();
	double EOSmuB();
	double EOSmuS();
	double EOSmuQ();

	double EOSp();
	double EOSs();
	double EOSe();
	double EOSw();
	double EOSs_terms_T(double Tin);
	//double EOScs2out();
	//double EOSwfz();
	//double EOScs2out();
	//double EOSwfz();
	double EOSA();
	double EOSdwds();
	double EOSs_out(double e_In);
	double EOSs_out(double e_In, double rhoB_In, double rhoS_In, double rhoQ_In);
	void EOSupdate_s(double s_In);
	void EOSupdate_s(double s_In, double rhoB_In, double rhoS_In, double rhoQ_In);

	double particle_T, particle_muB, particle_muS, particle_muQ;

};


#endif
