#ifndef _BBMG_H_
#define _BBMG_H_

#include<fstream>
#include<stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <cstdlib>
#include<iostream>
#include<string.h>
#include<dirent.h>
#include <sstream>
#include <vector>
#include "particle.h"
#include "LinkList.h"
#include <errno.h>
#include <sys/stat.h>
#include <string>
//#include <random>
#include <ctime>

using namespace std;

template <int D>
class BBMG {
private:
    struct field
    {
        int sph,on;
        double rho,rho0,T,v[2];
        double r[2],phi,line;
        int pid;
        double gam,vmag,vang;
    };

    double phi[15];
    double pre,pre2;
    double area;
    vector<double> rr;
    vector <field> ff;
    int z,a,c;
    double kappa,vjet;
    double Cg,Cq,q;
    double rho0tot; // total density, NOT just T>TD!!!
    double Rq[15],Rg[15];
    double flow(field &f);
    double efluc();
    double Pfg,Pfq;
    double gft(double p);
    double qft(double p);

public:
    BBMG<D>(LinkList<D> &linklist);
    void propogate(LinkList<D> &linklist);
    double kern(double r);
    void inter(LinkList<D> &linklist,field &f);
    void initial(LinkList<D> &linklist);
    double TD;
};


template <int D>
BBMG<D>::BBMG(LinkList<D> &linklist)
{
    TD=150;
    srand( time( NULL ) );
    q=1;
    Cg=3;
    Cq=4./3;
    //rr.resize(10);
    z=1; // path length dependence
    a=0; //factor for E
    c=2+z-a; //factor for T
    kappa=0.17;
    vjet=1;
    area=PI*pow(2.*linklist._h,2);
    rr.resize(linklist._n);

    pre=linklist.knorm;
    pre2=linklist.knorm2;

    for (int i=0; i<15; i++)
    {
        Rq[i]=0;
        Rg[i]=0;
        phi[i]=i*PI/7;
    }

    Pfg=10;
    Pfq=10;
}


//template <int D>
//double BBMG<D>::kappa(double temp)
//{
//    Te=270
//    Cr=2.134;
//    double kp=exp(-b*(T2-TD));
//    return (1-a)*Cr*kp;
//}


template <int D>
void BBMG<D>::initial(LinkList<D> &linklist)
{

    rho0tot=0;
    for (int i=0; i<linklist._n; i++)
    {
        double rsub=linklist._p[i].EOS.p()/linklist._p[i].EOS.T();
        rho0tot+=rsub;

        if (linklist._p[i].EOS.T()*197.3>TD)
        {
            field sub;
            sub.r[0]=linklist._p[i].r.x[0];
            sub.r[1]=linklist._p[i].r.x[1];
            sub.rho0=rsub;
            sub.sph=i;
            sub.line=0.5*pow(linklist.t0,z)*pow(sub.rho0,c)*linklist.dt; // only if initial flow=0

            for (int j=0; j<14; j++)
            {
                sub.phi=phi[j];
                sub.pid=j;
                sub.on=1;
                //sub.line=0.5*pow(linklist.t0,z)*pow(sub.rho,c)*flow(phi[j])*linklist.dt;
                ff.push_back(sub);
            }
        }
    }
}


template <int D>
double BBMG<D>::flow(field &f)
{
    return f.gam*(1-f.vmag*cos(f.phi-f.vang));
}


template <int D>
double BBMG<D>::gft(double p)
{
    return 2*p;
}


template <int D>
double BBMG<D>::qft(double p)
{
    return 2*p;
}


template <int D>
double BBMG<D>::efluc()
{
    int random_variable = std::rand()/RAND_MAX;
    double zeta=random_variable*(q+2.);

    return (1.+q)/pow(q+2,1+q)*pow(q+2.-zeta,q);
}


template <int D>
void BBMG<D>::propogate(LinkList<D> &linklist)
{
    double tau=linklist.t+linklist.t0;
    int stillon=0;
    int tot=ff.size();

    for (int i=0; i<tot; i++)
    {
        //propogate x,y position of jet
        ff[i].r[0]+=vjet*linklist.dt*cos(ff[i].phi);
        ff[i].r[1]+=vjet*linklist.dt*cos(ff[i].phi);

        inter(linklist,ff[i]);

        if ((ff[i].on==1)&&(ff[i].T>TD))
        {
            ff[i].line+=pow(tau,z)*pow(ff[i].rho,c)*flow(ff[i])*linklist.dt;
            stillon++;
        }
        else
        {
            ff[i].on=0;
            ff[i].line+=0.5*pow(tau,z)*pow(ff[i].rho,c)*flow(ff[i])*linklist.dt;
            ff[i].line*=kappa*efluc();

            double P0g=Pfg+Cg*pow(Pfg,1-a)*ff[i].line;
            double P0q=Pfq+Cq*pow(Pfq,1-a)*ff[i].line;

            int jj=ff[i].pid;

            Rq[jj]+=pow(P0g/Pfg,1+a)*ff[i].rho0*gft(P0g)/gft(Pfg);
            Rq[jj]+=pow(P0q/Pfq,1+a)*ff[i].rho0*qft(P0g)/qft(Pfg);
        }
    }

    if(stillon==0)
    {
        for (int j=0; j<14; j++)
        {
            Rq[j]/=rho0tot;
            Rg[j]/=rho0tot;
        }
    }


//    double maxx=0,maxy=0,minx=0,miny=0;
//    double rmax=0;
//    for (int i=0;i<linklist._n;i++)
//      {
//
//
//          rr[i]=sqrt(linklist._p[i].r.x[0]*linklist._p[i].r.x[0]+linklist._p[i].r.x[1]*linklist._p[i].r.x[1]);
//
//          if (linklist._p[i].EOS.T()*197.3>TD){
//              if(linklist._p[i].r.x[0]>maxx) maxx=linklist._p[i].r.x[0];
//              if(linklist._p[i].r.x[0]<minx) minx=linklist._p[i].r.x[0];
//              if(linklist._p[i].r.x[1]>maxy) maxy=linklist._p[i].r.x[1];
//              if(linklist._p[i].r.x[1]<miny) miny=linklist._p[i].r.x[1];
//              if(rr[i]>rmax) rmax=rr[i];
//          }
//      }
//
//
//

//    cout << rmax << " " << maxx << " " << maxy << " " << minx << " " << miny << endl;
}

template <int D>
void BBMG<D>::inter(LinkList<D> &linklist,field &f)
{

    double den=0,den2=0;
    for (int i=0; i<linklist._n; i++)
    {
        double dx=linklist._p[i].r.x[0]-f.r[0];
        double dy=linklist._p[i].r.x[1]-f.r[1];

        double rdiff=sqrt(dx*dx+dy*dy)/linklist._h;

        if (rdiff<2)
        {
            den++;
            den2+=linklist._p[i].sigmaweight;
            double kk=kern(rdiff);
            f.T+=linklist._p[i].EOS.T()*0.06*0.06*kk;
            f.rho+=(linklist._p[i].EOS.p()/linklist._p[i].EOS.T())*kk;
            f.v[0]+=linklist._p[i].v.x[0]*kk;
            f.v[1]+=linklist._p[i].v.x[1]*kk;

            cout << dx << " " << dy << " " << linklist._p[i].EOS.T()*197.3 << " " << linklist._p[i].v <<  endl;
        }
    }

    double fac=den/area;
    f.T*=197.3;
    f.rho/=fac;
    f.v[0]/=fac;
    f.v[1]/=fac;
    f.vmag=sqrt(pow(f.v[0],2)+pow(f.v[1],2));
    f.vang=atan2(f.v[1],f.v[0]);
    f.gam=1./sqrt( f.vmag*f.vmag + 1 );

//      cout << f.r[0] << " " << f.r[1] << " " << f.T << " " << f.v[0] << " " << f.v[1] << endl;
//      getchar();


}

template <int D>
double BBMG<D>::kern(double r)
{
    if(r>=2) return 0;

    if(r>=1) return pre2*pow((2.-r),3);

    return pre*(1. - 1.5*pow(r,2) + 0.75*pow(r,3));
}

#endif
