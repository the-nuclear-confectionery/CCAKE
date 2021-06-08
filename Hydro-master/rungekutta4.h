#ifndef ENTERIC_H_
#define ENTERIC_H_

#include <math.h>
#include "vector.h"
#include "matrix.h"
#include "LinkList.h"
#include "tables.h"
#include "particle.h"
#include "eos.h"

template <int D>
void adams(double dx,int & count,void (*derivatives)( LinkList<D> &linklist), LinkList<D> &linklist )
{
    int N=linklist.n();
    Vector<double,D> *u0,*r0;
    u0=new Vector<double,D> [N];
    r0=new Vector<double,D> [N];
    double E0,t0;

    linklist.rk2=1;
    if (count==0)
    {
        t0=linklist.t;
        for (int i=0; i<N; ++i) {
            u0[i]=linklist._p[i].u;
            r0[i]=linklist._p[i].r;
        }
        E0=linklist.Ez;


        //    first step
        (*derivatives)(linklist);
        for (int i=0; i<N; ++i) {

            linklist._p[i].u=u0[i]+0.5*dx*linklist._p[i].du_dt;
            linklist._p[i].r=r0[i]+0.5*dx*linklist._p[i].v;

            linklist._p[i].k1=linklist._p[i].du_dt;
            linklist._p[i].r1=linklist._p[i].v;
        }
        linklist.Ez=E0+0.5*dx*linklist.dEz;
        linklist.E1=linklist.dEz;
        linklist.t=t0+0.5*dx;


        //    second step
        (*derivatives)(linklist);

        for (int i=0; i<N; ++i) {
            linklist._p[i].u=u0[i]+    dx*linklist._p[i].du_dt;

            linklist._p[i].r=r0[i]+dx*linklist._p[i].v;

        }
        linklist.Ez=E0+dx*linklist.dEz;
        linklist.t=t0+dx;


        count=1;
    }
    else
    {
        t0=linklist.t;
        for (int i=0; i<N; ++i) {
            u0[i]=linklist._p[i].u;
            r0[i]=linklist._p[i].r;
        }
        E0=linklist.Ez;



        (*derivatives)(linklist);
        for (int i=0; i<N; ++i) {
            linklist._p[i].u=u0[i]+1.5*dx*linklist._p[i].du_dt-0.5*dx*linklist._p[i].k1;
            linklist._p[i].r=r0[i]+1.5*dx*linklist._p[i].v-0.5*dx*linklist._p[i].r1;

            linklist._p[i].k1=linklist._p[i].du_dt;
            linklist._p[i].r1=linklist._p[i].v;
        }
        linklist.Ez=E0+1.5*dx*linklist.dEz-0.5*dx*linklist.E1;
        linklist.t=t0+dx;


        linklist.E1=linklist.dEz;


    }

//    for (int i=0; i<linklist.n(); ++i) {
//    cout << dx << endl;
//    cout << linklist._p[i].u <<  "  " << u0[i] << " " << linklist.ux(linklist.t, linklist._p[i].r.x[0], linklist._p[i].r.x[1]) << " " << linklist.uy(linklist.t, linklist._p[i].r.x[0], linklist._p[i].r.x[1]) << " " << linklist._p[i].du_dt << endl;
//    getchar();
//    }

    delete [] u0;
    delete [] r0;

}

template <int D>
void adams2(double dx,int & count,void (*derivatives)( LinkList<D> &linklist), LinkList<D> &linklist )
{
    int N=linklist.n();
    Vector<double,D> *u0,*r0;
    u0=new Vector<double,D> [N];
    r0=new Vector<double,D> [N];
    double E0,t0;

    linklist.rk2=1;
    if (count==0)
    {
        t0=linklist.t;
        for (int i=0; i<N; ++i) {
            u0[i]=linklist._p[i].u;
            r0[i]=linklist._p[i].r;
        }
        E0=linklist.Ez;
        (*derivatives)(linklist);
        for (int i=0; i<N; ++i) {
            linklist._p[i].u=u0[i]+dx*linklist._p[i].du_dt;
            linklist._p[i].r=r0[i]+dx*linklist._p[i].v;

            linklist._p[i].k1=linklist._p[i].du_dt;
            linklist._p[i].r1=linklist._p[i].v;
        }
        linklist.Ez=E0+dx*linklist.dEz;
        linklist.E1=linklist.dEz;
        linklist.t=t0+dx;


        count=1;
    }
    else if (count==1)
    {
        t0=linklist.t;
        for (int i=0; i<N; ++i) {
            u0[i]=linklist._p[i].u;
            r0[i]=linklist._p[i].r;
        }
        E0=linklist.Ez;

        (*derivatives)(linklist);
        for (int i=0; i<N; ++i) {
            linklist._p[i].u=u0[i]+1.5*dx*linklist._p[i].du_dt-0.5*dx*linklist._p[i].k1;
            linklist._p[i].r=r0[i]+1.5*dx*linklist._p[i].v-0.5*dx*linklist._p[i].r1;

            linklist._p[i].k2=linklist._p[i].k1;
            linklist._p[i].r2=linklist._p[i].r1;

            linklist._p[i].k1=linklist._p[i].du_dt;
            linklist._p[i].r1=linklist._p[i].v;
        }
        linklist.Ez=E0+1.5*dx*linklist.dEz-0.5*dx*linklist.E1;
        linklist.t=t0+dx;

        linklist.E2=linklist.E1;
        linklist.E1=linklist.dEz;

        count=2;


    }
    else
    {
        t0=linklist.t;
        for (int i=0; i<N; ++i) {
            u0[i]=linklist._p[i].u;
            r0[i]=linklist._p[i].r;
        }
        E0=linklist.Ez;

        (*derivatives)(linklist);
        for (int i=0; i<N; ++i) {
            linklist._p[i].u=u0[i]+dx*(23./12.*linklist._p[i].du_dt-4./3.*linklist._p[i].k1+5/12.*linklist._p[i].k2);
            linklist._p[i].r=r0[i]+dx*(23./12.*linklist._p[i].v-4./3.*linklist._p[i].r1+5/12.*linklist._p[i].r2);

            linklist._p[i].k2=linklist._p[i].k1;
            linklist._p[i].r2=linklist._p[i].r1;

            linklist._p[i].k1=linklist._p[i].du_dt;
            linklist._p[i].r1=linklist._p[i].v;
        }
        linklist.Ez=E0+dx*(23./12.*linklist.dEz-4./3.*linklist.E1+5/12.*linklist.E2);
        linklist.t=t0+dx;

        linklist.E2=linklist.E1;
        linklist.E1=linklist.dEz;



    }

    delete [] u0;
    delete [] r0;

}


template <int D>
void rungeKutta2(double dx,void (*derivatives)( LinkList<D> &linklist), LinkList<D> &linklist )
{
    // creating arrays of vectors of the derivatives at each step
    int N=linklist.n();
    Vector<double,D> *u0,*r0;
    u0=new Vector<double,D> [N];
    r0=new Vector<double,D> [N];


    linklist.rk2=1;
    double t0=linklist.t;
    for (int i=0; i<N; ++i) {
        u0[i]=linklist._p[i].u;
        r0[i]=linklist._p[i].r;
    }
    double E0=linklist.Ez;

    //    first step
    (*derivatives)(linklist);
    for (int i=0; i<N; ++i) {
        linklist._p[i].u=u0[i]+0.5*dx*linklist._p[i].du_dt;

        linklist._p[i].r=r0[i]+0.5*dx*linklist._p[i].v;

    }
    linklist.Ez=E0+0.5*dx*linklist.dEz;
    linklist.t=t0+0.5*dx;

    //    second step
    (*derivatives)(linklist);
    for (int i=0; i<N; ++i) {
        linklist._p[i].u=u0[i]+    dx*linklist._p[i].du_dt;

        linklist._p[i].r=r0[i]+dx*linklist._p[i].v;

    }
    linklist.Ez=E0+dx*linklist.dEz;
    linklist.t=t0+dx;

    delete [] u0;
    delete [] r0;
}

template <int D>
void rungeKutta4(double dx,void (*derivatives)( LinkList<D> &linklist), LinkList<D> &linklist )
{
    // creating arrays of vectors of the derivatives at each step
    int N=linklist.n();
    Vector<double,D> *u0,*r0;
    u0=new Vector<double,D> [N];
    r0=new Vector<double,D> [N];
    double E0,E1,E2,E3,E4,t0;
//    double *etasig0;
//    etasig0=new double [N];
    //double Esub;

    linklist.rk2=1;
    t0=linklist.t;
    for (int i=0; i<N; ++i) {
        u0[i]=linklist._p[i].u;
        r0[i]=linklist._p[i].r;
//        etasig0[i]=linklist._p[i].eta_sigma;
    }
    E0=linklist.Ez;

    //    first step
    (*derivatives)(linklist);
    for (int i=0; i<N; ++i) {
        linklist._p[i].k1=dx*linklist._p[i].du_dt;
        linklist._p[i].u=u0[i]+0.5*linklist._p[i].k1;

        linklist._p[i].r1=dx*linklist._p[i].v;
        linklist._p[i].r=r0[i]+0.5*linklist._p[i].r1;

//        linklist._p[i].es1=dx*linklist._p[i].detasigma_dt;
//        linklist._p[i].eta_sigma=es0[i]+(1./2.0)*linklist._p[i].es1;

    }
    E1=dx*linklist.dEz;
    linklist.t=t0+0.5*dx;

    //    second step
    (*derivatives)(linklist);
    for (int i=0; i<N; ++i) {
        linklist._p[i].k2=dx*linklist._p[i].du_dt;
        linklist._p[i].u=u0[i]+0.5*linklist._p[i].k2;

        linklist._p[i].r2=dx*linklist._p[i].v;
        linklist._p[i].r=r0[i]+0.5*linklist._p[i].r2;

//        linklist._p[i].es2=dx*linklist._p[i].detasigma_dt;
//        linklist._p[i].eta_sigma=es0[i]+(1./2.0)*linklist._p[i].es2;
    }
    E2=dx*linklist.dEz;

    //    third step
    (*derivatives)(linklist);
    for (int i=0; i<N; ++i) {
        linklist._p[i].k3=dx*linklist._p[i].du_dt;
        linklist._p[i].u=u0[i]+linklist._p[i].k3;

        linklist._p[i].r3=dx*linklist._p[i].v;
        linklist._p[i].r=r0[i]+linklist._p[i].r3;

//        linklist._p[i].es3=dx*linklist._p[i].detasigma_dt;
//        linklist._p[i].eta_sigma=es0[i]+(1./2.0)*linklist._p[i].es3;
    }
    E3=dx*linklist.dEz;
    linklist.t=t0+dx;

    //    fourth step
    (*derivatives)(linklist);
    for (int i=0; i<N; ++i) {
        linklist._p[i].k4=dx*linklist._p[i].du_dt;

        linklist._p[i].r4=dx*linklist._p[i].v;

//        linklist._p[i].es4=dx*linklist._p[i].detasigma_dt;

        // sum the weighted steps into yf and return the final y values
        linklist._p[i].u=u0[i]+(1./6.0)*linklist._p[i].k1+(1./3.0)*linklist._p[i].k2+(1./3.0)*linklist._p[i].k3+(1./6.0)*linklist._p[i].k4;

        linklist._p[i].r=r0[i]+(1./6.0)*linklist._p[i].r1+(1./3.0)*linklist._p[i].r2+(1./3.0)*linklist._p[i].r3+(1./6.0)*linklist._p[i].r4;



//        linklist._p[i].eta_sigma=es0[i]+(1./6.0)*linklist._p[i].es1+(1./3.0)*linklist._p[i].es2+(1./3.0)*linklist._p[i].es3+(1./6.0)*linklist._p[i].es4;
    }
    E4=dx*linklist.dEz;
    linklist.Ez=E0+E1/6.+E2/3.+E3/3.+E4/6.;

    delete [] u0;
    delete [] r0;

}

/// viscous

template <int D>
void vadams(double dx,int & count,void (*derivatives)( LinkList<D> &linklist), LinkList<D> &linklist )
{
    int N=linklist.n();
    Vector<double,D> *u0,*r0;
    double *etasigma0,*Bulk0;
    u0=new Vector<double,D> [N];
    r0=new Vector<double,D> [N];
    etasigma0=new double [N];
    Bulk0=new double [N];
    double E0,t0;

    linklist.rk2=1;

    if (count==0)
    {
        t0=linklist.t;
        for (int i=0; i<N; ++i) {
            u0[i]=linklist._p[i].u;
            r0[i]=linklist._p[i].r;
            etasigma0[i]=linklist._p[i].eta_sigma;
            Bulk0[i]=linklist._p[i].Bulk;
        }
        E0=linklist.Ez;

        //    first step
        (*derivatives)(linklist);
        for (int i=0; i<N; ++i) {
            linklist._p[i].u=u0[i]+0.5*dx*linklist._p[i].du_dt;
            linklist._p[i].r=r0[i]+0.5*dx*linklist._p[i].v;
            linklist._p[i].eta_sigma=etasigma0[i]+0.5*dx*linklist._p[i].detasigma_dt;
            linklist._p[i].Bulk=Bulk0[i]+0.5*dx*linklist._p[i].dBulk_dt;

            linklist._p[i].k1=linklist._p[i].du_dt;
            linklist._p[i].r1=linklist._p[i].v;
            linklist._p[i].ets1=linklist._p[i].detasigma_dt;
            linklist._p[i].b1=linklist._p[i].dBulk_dt;
        }
        linklist.Ez=E0+0.5*dx*linklist.dEz;
        linklist.t=t0+0.5*dx;
        linklist.E1=linklist.dEz;

        //    second step
        (*derivatives)(linklist);
        for (int i=0; i<N; ++i) {
            linklist._p[i].u=u0[i]+    dx*linklist._p[i].du_dt;
            linklist._p[i].r=r0[i]+dx*linklist._p[i].v;
            linklist._p[i].eta_sigma=etasigma0[i]+dx*linklist._p[i].detasigma_dt;
            linklist._p[i].Bulk=Bulk0[i]+dx*linklist._p[i].dBulk_dt;

        }
        linklist.Ez=E0+0.5*dx*linklist.dEz+0.5*dx*linklist.dEz;
        linklist.t=t0+dx;


        count=1;
    }
    else
    {
        t0=linklist.t;
        for (int i=0; i<N; ++i) {
            u0[i]=linklist._p[i].u;
            r0[i]=linklist._p[i].r;
            etasigma0[i]=linklist._p[i].eta_sigma;
            Bulk0[i]=linklist._p[i].Bulk;
        }
        E0=linklist.Ez;

        (*derivatives)(linklist);
        for (int i=0; i<N; ++i) {
            linklist._p[i].u=u0[i]+1.5*dx*linklist._p[i].du_dt-0.5*dx*linklist._p[i].k1;
            linklist._p[i].r=r0[i]+1.5*dx*linklist._p[i].v-0.5*dx*linklist._p[i].r1;
            linklist._p[i].eta_sigma=etasigma0[i]+1.5*dx*linklist._p[i].detasigma_dt-0.5*dx*linklist._p[i].ets1;
            linklist._p[i].Bulk=Bulk0[i]+1.5*dx*linklist._p[i].dBulk_dt-0.5*dx*linklist._p[i].b1;

            linklist._p[i].k1=linklist._p[i].du_dt;
            linklist._p[i].r1=linklist._p[i].v;
            linklist._p[i].ets1=linklist._p[i].detasigma_dt;
            linklist._p[i].b1=linklist._p[i].dBulk_dt;
        }
        linklist.Ez=E0+1.5*dx*linklist.dEz-0.5*dx*linklist.E1;
        linklist.t=t0+dx;

        linklist.E1=linklist.dEz;


    }

    delete [] u0;
    delete [] r0;
    delete [] etasigma0;
    delete [] Bulk0;
}




template <int D>
void vrungeKutta2(double dx,void (*derivatives)( LinkList<D> &linklist), LinkList<D> &linklist )
{
    // creating arrays of vectors of the derivatives at each step
    int N=linklist.n();
    Vector<double,D> *u0,*r0;
    double *etasigma0,*Bulk0;
    u0=new Vector<double,D> [N];
    r0=new Vector<double,D> [N];
    etasigma0=new double [N];
    Bulk0=new double [N];
    double E0,t0;

    linklist.rk2=1;
    t0=linklist.t;
    for (int i=0; i<N; ++i) {
        u0[i]=linklist._p[i].u;
        r0[i]=linklist._p[i].r;
        etasigma0[i]=linklist._p[i].eta_sigma;
        Bulk0[i]=linklist._p[i].Bulk;
    }
    E0=linklist.Ez;

    //    first step
    (*derivatives)(linklist);
    for (int i=0; i<N; ++i) {
        linklist._p[i].u=u0[i]+0.5*dx*linklist._p[i].du_dt;
        linklist._p[i].r=r0[i]+0.5*dx*linklist._p[i].v;
        linklist._p[i].eta_sigma=etasigma0[i]+0.5*dx*linklist._p[i].detasigma_dt;
        linklist._p[i].Bulk=Bulk0[i]+0.5*dx*linklist._p[i].dBulk_dt;
        //if (linklist._p[i].eta_sigma<0) linklist._p[i].eta_sigma=0;

        //if (i==14335) cout <<linklist._p[i].eta_sigma << " " << linklist._p[i].u << " " <<  linklist._p[i].Bulk << endl;
    }
    linklist.Ez=E0+0.5*dx*linklist.dEz;
    linklist.t=t0+0.5*dx;

    //    second step
    (*derivatives)(linklist);
    for (int i=0; i<N; ++i) {
        linklist._p[i].u=u0[i]+    dx*linklist._p[i].du_dt;
        linklist._p[i].r=r0[i]+dx*linklist._p[i].v;
        linklist._p[i].eta_sigma=etasigma0[i]+dx*linklist._p[i].detasigma_dt;
        linklist._p[i].Bulk=Bulk0[i]+dx*linklist._p[i].dBulk_dt;
        //if (linklist._p[i].eta_sigma<0) linklist._p[i].eta_sigma=0;

        //if (i==14335) cout <<linklist._p[i].eta_sigma << " " << linklist._p[i].u << " " << linklist._p[i].Bulk << endl;

    }
    linklist.Ez=E0+dx*linklist.dEz;
    linklist.t=t0+dx;

    delete [] u0;
    delete [] r0;
    delete [] etasigma0;
    delete [] Bulk0;
}

template <int D>
void vrungeKutta4(double dx,void (*derivatives)( LinkList<D> &linklist), LinkList<D> &linklist )
{
    // creating arrays of vectors of the derivatives at each step
    int N=linklist.n();
    Vector<double,D> *u0,*r0;
    double *etasigma0,*Bulk0;
    u0=new Vector<double,D> [N];
    r0=new Vector<double,D> [N];
    etasigma0=new double [N];
    Bulk0=new double [N];
    double E0,E1,E2,E3,E4,t0;

    linklist.rk2=1;
    t0=linklist.t;
    for (int i=0; i<N; ++i) {
        u0[i]=linklist._p[i].u;
        r0[i]=linklist._p[i].r;
        etasigma0[i]=linklist._p[i].eta_sigma;
        Bulk0[i]=linklist._p[i].Bulk;
    }
    E0=linklist.Ez;

    //    first step

    (*derivatives)(linklist);

    for (int i=0; i<N; ++i) {
        linklist._p[i].k1=dx*linklist._p[i].du_dt;
        linklist._p[i].u=u0[i]+0.5*linklist._p[i].k1;

        linklist._p[i].r1=dx*linklist._p[i].v;
        linklist._p[i].r=r0[i]+0.5*linklist._p[i].r1;

        linklist._p[i].ets1=dx*linklist._p[i].detasigma_dt;
        linklist._p[i].eta_sigma=etasigma0[i]+0.5*linklist._p[i].ets1;

        linklist._p[i].b1=dx*linklist._p[i].dBulk_dt;
        linklist._p[i].Bulk=Bulk0[i]+0.5*linklist._p[i].b1;


    }
    E1=dx*linklist.dEz;
    linklist.t=t0+0.5*dx;


    //    second step
    (*derivatives)(linklist);
    for (int i=0; i<N; ++i) {
        linklist._p[i].k2=dx*linklist._p[i].du_dt;
        linklist._p[i].u=u0[i]+0.5*linklist._p[i].k2;

        linklist._p[i].r2=dx*linklist._p[i].v;
        linklist._p[i].r=r0[i]+0.5*linklist._p[i].r2;

        linklist._p[i].ets2=dx*linklist._p[i].detasigma_dt;
        linklist._p[i].eta_sigma=etasigma0[i]+0.5*linklist._p[i].ets2;

        linklist._p[i].b2=dx*linklist._p[i].dBulk_dt;
        linklist._p[i].Bulk=Bulk0[i]+0.5*linklist._p[i].b2;

    }
    E2=dx*linklist.dEz;

    //    third step
    (*derivatives)(linklist);
    for (int i=0; i<N; ++i) {
        linklist._p[i].k3=dx*linklist._p[i].du_dt;
        linklist._p[i].u=u0[i]+linklist._p[i].k3;

        linklist._p[i].r3=dx*linklist._p[i].v;
        linklist._p[i].r=r0[i]+linklist._p[i].r3;

        linklist._p[i].ets3=dx*linklist._p[i].detasigma_dt;
        linklist._p[i].eta_sigma=etasigma0[i]+linklist._p[i].ets3;

        linklist._p[i].b3=dx*linklist._p[i].dBulk_dt;
        linklist._p[i].Bulk=Bulk0[i]+linklist._p[i].b3;

    }
    E3=dx*linklist.dEz;
    linklist.t=t0+dx;

    //    fourth step
    (*derivatives)(linklist);
    for (int i=0; i<N; ++i) {
        linklist._p[i].k4=dx*linklist._p[i].du_dt;

        linklist._p[i].r4=dx*linklist._p[i].v;

        linklist._p[i].ets4=dx*linklist._p[i].detasigma_dt;

        linklist._p[i].b4=dx*linklist._p[i].dBulk_dt;



        // sum the weighted steps into yf and return the final y values
        linklist._p[i].u=u0[i]+(1./6.0)*linklist._p[i].k1+(1./3.0)*linklist._p[i].k2+(1./3.0)*linklist._p[i].k3+(1./6.0)*linklist._p[i].k4;

        linklist._p[i].r=r0[i]+(1./6.0)*linklist._p[i].r1+(1./3.0)*linklist._p[i].r2+(1./3.0)*linklist._p[i].r3+(1./6.0)*linklist._p[i].r4;


        linklist._p[i].eta_sigma=etasigma0[i]+(1./6.0)*linklist._p[i].ets1+(1./3.0)*linklist._p[i].ets2+(1./3.0)*linklist._p[i].ets3+(1./6.0)*linklist._p[i].ets4;

        linklist._p[i].Bulk=Bulk0[i]+(1./6.0)*linklist._p[i].b1+(1./3.0)*linklist._p[i].b2+(1./3.0)*linklist._p[i].b3+(1./6.0)*linklist._p[i].b4;


    }
    E4=dx*linklist.dEz;
    linklist.Ez=E0+E1/6.+E2/3.+E3/3.+E4/6.;

    delete [] u0;
    delete [] r0;
    delete [] etasigma0;
    delete [] Bulk0;
}

// SHEAR VISCOSITY

template <int D>
void svadams(double dx,int & count,void (*derivatives)( LinkList<D> &linklist), LinkList<D> &linklist )
{
    int N=linklist.n();
    Vector<double,D> *u0,*r0;
    double *etasigma0,*Bulk0;
    Matrix<double,D,D> *shv0;
    u0=new Vector<double,D> [N];
    r0=new Vector<double,D> [N];
    etasigma0=new double [N];
    Bulk0=new double [N];
    shv0=new Matrix<double,D,D> [N];
    double E0,t0;

    linklist.rk2=1;

    if (count==0)
    {
        t0=linklist.t;
        for (int i=0; i<N; ++i) {
            u0[i]=linklist._p[i].u;
            r0[i]=linklist._p[i].r;
            etasigma0[i]=linklist._p[i].eta_sigma;
            Bulk0[i]=linklist._p[i].Bulk;
            mini(shv0[i],linklist._p[i].shv);
        }
        E0=linklist.Ez;

        //    first step
        (*derivatives)(linklist);
        for (int i=0; i<N; ++i) {
            linklist._p[i].u=u0[i]+0.5*dx*linklist._p[i].du_dt;
            linklist._p[i].r=r0[i]+0.5*dx*linklist._p[i].v;
            linklist._p[i].eta_sigma=etasigma0[i]+0.5*dx*linklist._p[i].detasigma_dt;
            linklist._p[i].Bulk=Bulk0[i]+0.5*dx*linklist._p[i].dBulk_dt;
            tmini(linklist._p[i].shv,shv0[i]+0.5*dx*linklist._p[i].dshv_dt);

            linklist._p[i].k1=linklist._p[i].du_dt;
            linklist._p[i].r1=linklist._p[i].v;
            linklist._p[i].ets1=linklist._p[i].detasigma_dt;
            linklist._p[i].b1=linklist._p[i].dBulk_dt;
            linklist._p[i].shv1=linklist._p[i].dshv_dt;
        }
        linklist.Ez=E0+0.5*dx*linklist.dEz;
        linklist.t=t0+0.5*dx;
        linklist.E1=linklist.dEz;

        //    second step
        (*derivatives)(linklist);
        for (int i=0; i<N; ++i) {
            linklist._p[i].u=u0[i]+    dx*linklist._p[i].du_dt;
            linklist._p[i].r=r0[i]+dx*linklist._p[i].v;
            linklist._p[i].eta_sigma=etasigma0[i]+dx*linklist._p[i].detasigma_dt;
            linklist._p[i].Bulk=Bulk0[i]+dx*linklist._p[i].dBulk_dt;
            tmini(linklist._p[i].shv,shv0[i]+dx*linklist._p[i].dshv_dt);

        }
        linklist.Ez=E0+0.5*dx*linklist.dEz+0.5*dx*linklist.dEz;
        linklist.t=t0+dx;


        count=1;
    }
    else
    {
        t0=linklist.t;
        for (int i=0; i<N; ++i) {
            u0[i]=linklist._p[i].u;
            r0[i]=linklist._p[i].r;
            etasigma0[i]=linklist._p[i].eta_sigma;
            Bulk0[i]=linklist._p[i].Bulk;
            shv0[i]=linklist._p[i].shv;
        }
        E0=linklist.Ez;

        (*derivatives)(linklist);
        for (int i=0; i<N; ++i) {
            linklist._p[i].u=u0[i]+1.5*dx*linklist._p[i].du_dt-0.5*dx*linklist._p[i].k1;
            linklist._p[i].r=r0[i]+1.5*dx*linklist._p[i].v-0.5*dx*linklist._p[i].r1;
            linklist._p[i].eta_sigma=etasigma0[i]+1.5*dx*linklist._p[i].detasigma_dt-0.5*dx*linklist._p[i].ets1;
            linklist._p[i].Bulk=Bulk0[i]+1.5*dx*linklist._p[i].dBulk_dt-0.5*dx*linklist._p[i].b1;
            tmini(linklist._p[i].shv,shv0[i]+1.5*dx*linklist._p[i].dshv_dt-0.5*dx*linklist._p[i].shv1);

            linklist._p[i].k1=linklist._p[i].du_dt;
            linklist._p[i].r1=linklist._p[i].v;
            linklist._p[i].ets1=linklist._p[i].detasigma_dt;
            linklist._p[i].b1=linklist._p[i].dBulk_dt;
            linklist._p[i].shv1=linklist._p[i].dshv_dt;
        }
        linklist.Ez=E0+1.5*dx*linklist.dEz-0.5*dx*linklist.E1;
        linklist.t=t0+dx;

        linklist.E1=linklist.dEz;


    }

    delete [] u0;
    delete [] r0;
    delete [] etasigma0;
    delete [] Bulk0;
    delete [] shv0;
}




template <int D>
void svrungeKutta2(double dx,void (*derivatives)( LinkList<D> &linklist), LinkList<D> &linklist )
{
    // creating arrays of vectors of the derivatives at each step
    int N=linklist.n();
    Vector<double,D> *u0,*r0;
    double *etasigma0,*Bulk0;
    u0=new Vector<double,D> [N];
    r0=new Vector<double,D> [N];
    etasigma0=new double [N];
    Bulk0=new double [N];
    Matrix <double,D,D> *shv0;
    shv0=new Matrix <double,D,D> [N];
    double E0,t0;




    linklist.rk2=1;
    t0=linklist.t;



    for (int i=0; i<N; ++i) {
        u0[i]=linklist._p[i].u;
        r0[i]=linklist._p[i].r;
        etasigma0[i]=linklist._p[i].eta_sigma;
        Bulk0[i]=linklist._p[i].Bulk;
        mini(shv0[i],linklist._p[i].shv);
    }

    E0=linklist.Ez;


    //    first step
    (*derivatives)(linklist);
    for (int i=0; i<N; ++i) {
        linklist._p[i].u=u0[i]+0.5*dx*linklist._p[i].du_dt;
        linklist._p[i].r=r0[i]+0.5*dx*linklist._p[i].v;
        linklist._p[i].eta_sigma=etasigma0[i]+0.5*dx*linklist._p[i].detasigma_dt;
        linklist._p[i].Bulk=Bulk0[i]+0.5*dx*linklist._p[i].dBulk_dt;
        tmini(linklist._p[i].shv,shv0[i]+0.5*dx*linklist._p[i].dshv_dt);


        //if (i==14335) cout <<linklist._p[i].shv << " " << shv0[i]+0.5*dx*linklist._p[i].dshv_dt << endl;


        //if (linklist._p[i].eta_sigma<0) linklist._p[i].eta_sigma=0;
    }
    linklist.Ez=E0+0.5*dx*linklist.dEz;
    linklist.t=t0+0.5*dx;

    //    second step
    (*derivatives)(linklist);
    for (int i=0; i<N; ++i) {
        linklist._p[i].u=u0[i]+    dx*linklist._p[i].du_dt;
        linklist._p[i].r=r0[i]+dx*linklist._p[i].v;
        linklist._p[i].eta_sigma=etasigma0[i]+dx*linklist._p[i].detasigma_dt;
        linklist._p[i].Bulk=Bulk0[i]+dx*linklist._p[i].dBulk_dt;
        tmini(linklist._p[i].shv,shv0[i]+dx*linklist._p[i].dshv_dt);
        //if (linklist._p[i].eta_sigma<0) linklist._p[i].eta_sigma=0;

        //if (i==14335) cout <<linklist._p[i].eta_sigma << " " << linklist._p[i].u << " " <<  linklist._p[i].Bulk << endl;
    }
    linklist.Ez=E0+dx*linklist.dEz;
    linklist.t=t0+dx;

    delete [] u0;
    delete [] r0;
    delete [] etasigma0;
    delete [] Bulk0;
    delete [] shv0;
}

template <int D>
void svrungeKutta4(double dx,void (*derivatives)( LinkList<D> &linklist), LinkList<D> &linklist )
{
    // creating arrays of vectors of the derivatives at each step
    int N=linklist.n();
    Vector<double,D> *u0,*r0;
    double *etasigma0,*Bulk0;
    u0=new Vector<double,D> [N];
    r0=new Vector<double,D> [N];
    etasigma0=new double [N];
    Bulk0=new double [N];
    Matrix<double,D,D> *shv0;
    shv0=new Matrix<double,D,D> [N];
    double E0,E1,E2,E3,E4,t0;

    linklist.rk2=1;
    t0=linklist.t;
    for (int i=0; i<N; ++i) {
        u0[i]=linklist._p[i].u;
        r0[i]=linklist._p[i].r;
        etasigma0[i]=linklist._p[i].eta_sigma;
        Bulk0[i]=linklist._p[i].Bulk;
        mini(shv0[i],linklist._p[i].shv);
    }
    E0=linklist.Ez;

    //    first step

    (*derivatives)(linklist);

    for (int i=0; i<N; ++i) {
        linklist._p[i].k1=dx*linklist._p[i].du_dt;
        linklist._p[i].u=u0[i]+0.5*linklist._p[i].k1;

        linklist._p[i].r1=dx*linklist._p[i].v;
        linklist._p[i].r=r0[i]+0.5*linklist._p[i].r1;

        linklist._p[i].ets1=dx*linklist._p[i].detasigma_dt;
        linklist._p[i].eta_sigma=etasigma0[i]+0.5*linklist._p[i].ets1;

        linklist._p[i].b1=dx*linklist._p[i].dBulk_dt;
        linklist._p[i].Bulk=Bulk0[i]+0.5*linklist._p[i].b1;

        linklist._p[i].shv1=dx*linklist._p[i].dshv_dt;
        tmini(linklist._p[i].shv,shv0[i]+0.5*linklist._p[i].shv1);


    }
    E1=dx*linklist.dEz;
    linklist.t=t0+0.5*dx;


    //    second step
    (*derivatives)(linklist);
    for (int i=0; i<N; ++i) {
        linklist._p[i].k2=dx*linklist._p[i].du_dt;
        linklist._p[i].u=u0[i]+0.5*linklist._p[i].k2;

        linklist._p[i].r2=dx*linklist._p[i].v;
        linklist._p[i].r=r0[i]+0.5*linklist._p[i].r2;

        linklist._p[i].ets2=dx*linklist._p[i].detasigma_dt;
        linklist._p[i].eta_sigma=etasigma0[i]+0.5*linklist._p[i].ets2;

        linklist._p[i].b2=dx*linklist._p[i].dBulk_dt;
        linklist._p[i].Bulk=Bulk0[i]+0.5*linklist._p[i].b2;

        linklist._p[i].shv2=dx*linklist._p[i].dshv_dt;
        tmini(linklist._p[i].shv,shv0[i]+0.5*linklist._p[i].shv2);

    }
    E2=dx*linklist.dEz;

    //    third step
    (*derivatives)(linklist);
    for (int i=0; i<N; ++i) {
        linklist._p[i].k3=dx*linklist._p[i].du_dt;
        linklist._p[i].u=u0[i]+linklist._p[i].k3;

        linklist._p[i].r3=dx*linklist._p[i].v;
        linklist._p[i].r=r0[i]+linklist._p[i].r3;

        linklist._p[i].ets3=dx*linklist._p[i].detasigma_dt;
        linklist._p[i].eta_sigma=etasigma0[i]+linklist._p[i].ets3;

        linklist._p[i].b3=dx*linklist._p[i].dBulk_dt;
        linklist._p[i].Bulk=Bulk0[i]+linklist._p[i].b3;

        linklist._p[i].shv3=dx*linklist._p[i].dshv_dt;
        tmini(linklist._p[i].shv,shv0[i]+0.5*linklist._p[i].shv3);

    }
    E3=dx*linklist.dEz;
    linklist.t=t0+dx;

    //    fourth step
    (*derivatives)(linklist);
    for (int i=0; i<N; ++i) {
        linklist._p[i].k4=dx*linklist._p[i].du_dt;

        linklist._p[i].r4=dx*linklist._p[i].v;

        linklist._p[i].ets4=dx*linklist._p[i].detasigma_dt;

        linklist._p[i].b4=dx*linklist._p[i].dBulk_dt;

        linklist._p[i].shv4=dx*linklist._p[i].dshv_dt;


        // sum the weighted steps into yf and return the final y values
        linklist._p[i].u=u0[i]+(1./6.0)*linklist._p[i].k1+(1./3.0)*linklist._p[i].k2+(1./3.0)*linklist._p[i].k3+(1./6.0)*linklist._p[i].k4;

        linklist._p[i].r=r0[i]+(1./6.0)*linklist._p[i].r1+(1./3.0)*linklist._p[i].r2+(1./3.0)*linklist._p[i].r3+(1./6.0)*linklist._p[i].r4;


        linklist._p[i].eta_sigma=etasigma0[i]+(1./6.0)*linklist._p[i].ets1+(1./3.0)*linklist._p[i].ets2+(1./3.0)*linklist._p[i].ets3+(1./6.0)*linklist._p[i].ets4;

        linklist._p[i].Bulk=Bulk0[i]+(1./6.0)*linklist._p[i].b1+(1./3.0)*linklist._p[i].b2+(1./3.0)*linklist._p[i].b3+(1./6.0)*linklist._p[i].b4;

        tmini(linklist._p[i].shv,shv0[i]+(1./6.0)*linklist._p[i].shv1+(1./3.0)*linklist._p[i].shv2+(1./3.0)*linklist._p[i].shv3+(1./6.0)*linklist._p[i].shv4);
    }
    E4=dx*linklist.dEz;
    linklist.Ez=E0+E1/6.+E2/3.+E3/3.+E4/6.;

    delete [] u0;
    delete [] r0;
    delete [] etasigma0;
    delete [] Bulk0;
    delete [] shv0;
}

template <int D>
void bsqrungeKutta2(double dx,void (*derivatives)( LinkList<D> &linklist), LinkList<D> &linklist )
{
cout << "Made it to " << __PRETTY_FUNCTION__ << "::" << __LINE__ << "!" << endl;

    // creating arrays of vectors of the derivatives at each step
    int N=linklist.n();
    Vector<double,D> *u0,*r0;
    double *etasigma0,*Bulk0;
	//double *rhoB0,*rhoS0,*rhoQ0;
    u0=new Vector<double,D> [N];
    r0=new Vector<double,D> [N];
    etasigma0=new double [N];
    Bulk0=new double [N];
    //rhoB0=new double [N];
    //rhoS0=new double [N];
    //rhoQ0=new double [N];

    Matrix <double,D,D> *shv0;
    shv0=new Matrix <double,D,D> [N];
    double E0,t0;




    linklist.rk2=1;
    t0=linklist.t;



    for (int i=0; i<N; ++i) {
        u0[i]=linklist._p[i].u;
        r0[i]=linklist._p[i].r;
        etasigma0[i]=linklist._p[i].eta_sigma;
        Bulk0[i]=linklist._p[i].Bulk;
        //rhoB0[i]=linklist._p[i].rhoB;
        //rhoS0[i]=linklist._p[i].rhoS;
        //rhoQ0[i]=linklist._p[i].rhoQ;
        mini(shv0[i],linklist._p[i].shv);
    }

    E0=linklist.Ez;


    //    first step
    (*derivatives)(linklist);
    for (int i=0; i<N; ++i) {
        linklist._p[i].u=u0[i]+0.5*dx*linklist._p[i].du_dt;
        linklist._p[i].r=r0[i]+0.5*dx*linklist._p[i].v;
        linklist._p[i].eta_sigma=etasigma0[i]+0.5*dx*linklist._p[i].detasigma_dt;
        linklist._p[i].Bulk=Bulk0[i]+0.5*dx*linklist._p[i].dBulk_dt;
        //linklist._p[i].rhoB=rhoB0[i]+0.5*dx*linklist._p[i].drhoB_dt;
        //linklist._p[i].rhoS=rhoS0[i]+0.5*dx*linklist._p[i].drhoS_dt;
        //linklist._p[i].rhoQ=rhoQ0[i]+0.5*dx*linklist._p[i].drhoQ_dt;
        tmini(linklist._p[i].shv,shv0[i]+0.5*dx*linklist._p[i].dshv_dt);


        //if (i==14335) cout <<linklist._p[i].shv << " " << shv0[i]+0.5*dx*linklist._p[i].dshv_dt << endl;


        //if (linklist._p[i].eta_sigma<0) linklist._p[i].eta_sigma=0;
    }
    linklist.Ez=E0+0.5*dx*linklist.dEz;
    linklist.t=t0+0.5*dx;

    //    second step
    (*derivatives)(linklist);
    for (int i=0; i<N; ++i) {
        linklist._p[i].u=u0[i]+    dx*linklist._p[i].du_dt;
        linklist._p[i].r=r0[i]+dx*linklist._p[i].v;
        linklist._p[i].eta_sigma=etasigma0[i]+dx*linklist._p[i].detasigma_dt;
std::cout << "Check eta_sigma: " << i << "   "
			<< etasigma0[i] << "   "
			<< dx << "   "
			<< linklist._p[i].detasigma_dt << std::endl;
        //linklist._p[i].rhoB=rhoB0[i]+dx*linklist._p[i].drhoB_dt;
        //linklist._p[i].rhoS=rhoS0[i]+dx*linklist._p[i].drhoS_dt;
        //linklist._p[i].rhoQ=rhoQ0[i]+dx*linklist._p[i].drhoQ_dt;
        tmini(linklist._p[i].shv,shv0[i]+dx*linklist._p[i].dshv_dt);
        //if (linklist._p[i].eta_sigma<0) linklist._p[i].eta_sigma=0;

        //if (i==14335) cout <<linklist._p[i].eta_sigma << " " << linklist._p[i].u << " " <<  linklist._p[i].Bulk << endl;
    }
    linklist.Ez=E0+dx*linklist.dEz;
    linklist.t=t0+dx;

    delete [] u0;
    delete [] r0;
    delete [] etasigma0;
    delete [] Bulk0;
    //delete [] rhoB0;
    //delete [] rhoS0;
    //delete [] rhoQ0;
    delete [] shv0;
}


#endif
