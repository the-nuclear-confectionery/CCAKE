
#ifndef _LINKLIST_H_
#define _LINKLIST_H_

#include "vector.h"
#include "matrix.h"
#include "particle.h"
#include "mathdef.h"
#include "random.h"
#include "eos.h"
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <vector>

extern char ifolder [];

template <int D>
class LinkList {
private:


    Vector<double,D> min;
    Vector<double,D> max;
    int range; //range is number of boxes from left ot right extra
    int Size;
    static constexpr double e0=1.;
    static constexpr double q=1.;
    double step;
    static constexpr double dTemp=0.00001;



    Vector<double,D> uni;


    void qmflow();



    double gradPressure_weight(int a, int b) {

        double _alpha_q = 1.;
        double _v_signal_q = sqrt(1./3.);

        double innerp = inner(_p[a].r -_p[b].r,_p[a].qmom -_p[b].qmom);
        double innerr = inner(_p[a].r -_p[b].r,_p[a].r-_p[b].r);
        innerp=2.*_alpha_q*_v_signal_q/(_p[a].sigma/_p[a].gamma+_p[b].sigma/_p[b].gamma)/sqrt(innerr)*innerp;

        if(innerp>0. || a==b) innerp=0.;


        return _p[b].sigmaweight*_p[a].sigma*( _p[b].EOSp()/(_p[b].sigma*_p[b].sigma) + _p[a].EOSp()/(_p[a].sigma*_p[a].sigma) -innerp);
    }
//    double gradBulk_weight(int a, int b) {    return _p[b].sigmaweight*_p[a].sigma*( (_p[b].Bulk_N+_p[b].Bulk)/_p[b].sigma/_p[b].gamma/t + (_p[a].Bulk_N+_p[a].Bulk)/_p[a].sigma/_p[a].gamma/t );}

    Vector<double,D> gradKernel (Vector<double,D> a);

    double linint(double x, double x1, double x2, double y1, double y2);
    Vector<double,D> linint(double x, double x1, double x2, Vector<double,D> y1, Vector<double,D> y2);
    double linfac(double fac, double y1, double y2);
    Vector<double,D> linfac(double fac, Vector<double,D> y1, Vector<double,D> y2);
    Matrix<double,D+1,D+1> linfac(double fac, Matrix<double,D+1,D+1> y1, Matrix<double,D+1,D+1> y2);
public:
    double _h;
    int _n;
    int start,end,fnum;
    vector<double> sFO,Tfluc; //entropy at freezeout
    int qmf; //if==1 quantum mechanicanical corrections to the flow or added, if==0 no corrections are included

    double kernel (Vector<double,D> kin);

    double knorm,knorm2,kgrad,kgrad2;
    int number_part;
    double t0;
    static constexpr double tend=50.02;
    double t,dt;
    double factor;
    int frzc;
    double tau,taup, taupp;
    int rk2;
    double gd2;

    int gtyp;
    //double van,avgT,avgvt;

    int cfon;
    vector<int> list;
    int *lead;
    int *link;
    int cf;
    int visc; //visc=0 for ideal, visc=1 for bulk, visc=2 for shear, visc=3 for bulk+shear
    double efcheck,sfcheck;
    int steps;
    //double *ep1,*ep2,*ec1,*ec2,*tsave,*pan,*san;

    int first;
    int average;
    int lowT;

    double *divTtemp,*gsub,*bulksub,*swsub,*shear33sub,*tlist;
    double avgetasig; // possibly not needed?
    Matrix<double,D+1,D+1> *shearsub;
    Vector<double,D> *divT,*rsub;
    Vector<double,D> *uout;
    double wfz,cs2;

    Vector<int,D> size;

    Vector<int,D> *dael;

    Particle<D> *_p;

    double dEz,E,Ez,E0,Etot,S,S0,Eloss,Esubb;
	double Btotal, Stotal, Qtotal;
	double Btotal0, Stotal0, Qtotal0;

    double E1,E2;

    LinkList<D>();

    int etaconst;
    double bvf, svf,zwidth,sTc,zTc;

    void setup(double it0, int ntot,double h,Particle<D> *_p,double dtsave, int & numpart);
    void setupnext(int ntot, Particle<D> *_pin, int & numpart);
    void vfreezeout(int curfrz);
    void svfreezeout(int curfrz);
    void bsqsvfreezeout(int curfrz);

    string eos_s,eos_p;
    string *filenames;
    string ebe_folder;

    int fcount,cevent;
    string eost;

    ~LinkList<D>();

    void initiate();

    void guess();
    void guess2();

    int triToSum(Vector<int,D> dael, Vector<int,D> size);

    void destroy();

    void optimization(int a);
    void optimization2(int a);
    void voptimization(int a);
    void voptimization2(int a,double tin);
    void svoptimization2(int a,double tin,int & count);
    void bsqsvoptimization(int a);
    void bsqsvoptimization2(int a,double tin,int & count);
    void conservation_entropy();
    void bsqconservation_entropy();
    void conservation();
    void conservation_E();
    void conservation_Ez();
    void vconservation();
    void vconservation_E();
    void vconservation_Ez();
    void svconservation();
    void svconservation_E();
    void svconservation_Ez();
    void bsqsvconservation();
    void bsqsvconservation_E();
    void bsqsvconservation_Ez();
	void conservation_BSQ();	// function to check conservation of B, S, and Q
    int n() {
        return _n;
    }
    void freezeout(int curfrz);
    //void vfreezeout(int curfrz);
    //void svfreezeout(int curfrz);
    //void bsqsvfreezeout(int curfrz);
    void interpolate(int curfrz);
    void vinterpolate(int curfrz);
    void svinterpolate(int curfrz);
    void bsqsvinterpolate(int curfrz);
    //void normals(int curfrz);
    void freezeset();
    void bsqsvfreezeset();
    void smoothopt(int a);
    void optint(int a, double & T00,  double & Tx0);


    void gubser(double h);
    void gubsershear( double h);
    void gubsershearbsq( double h);
    double eanal2(double t, double x, double y);
    double g(double t, double x, double y);
    double ux(double t, double x, double y);
    double uy(double t, double x, double y);
    void updateIC();
//    void p_anis(int &in);
//    void vp_anis(int &in);
//    void s_anis(int in);
//    void averages();
//    void v_anis();
    void setv(string vtype);
    void etas_set();
    void sv_set();
    void sv_setb();
    void bsqsv_set();
    void bsqsv_setb();
    void endEV();
//    void printsw();
    void setshear();
    void prints(); //possibly not needed?

};

#endif
