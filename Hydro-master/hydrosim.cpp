#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <sstream>

using namespace std;

//#include SPH.h
#include "vector.h"
#include "tables.h"
#include "particle.h"
#include "rungekutta4.h"
#include "hydrosim.h"
#include "eos.h"
#include "output.h"
#include "bbmg.h"


// here the hydro simulation is actually called within the various Simulation subroutines.  First the LinkList is initiated, then the SPH method is used, next we check if any particles have reached freezeout and finally we run Runge Kutta.  The hydro  equations of motion used within Runge Kutta are sent to the function and can be stored in any function you like.





void Simulation(double dt,LinkList<2> &linklist)
{
    cout << "Ready to start hydrodynamics\n";


    // int cc=0; // uncomment for adams method


    linklist.frzc=0;

    linklist.cf=0;


    Output<2> out(linklist); // sets up output

    linklist.t=linklist.t0;

    linklist.Ez=0;

    if ((linklist.qmf==1)||(linklist.qmf==3)) {
        out.eprofile(linklist);
        cout << "printed first timestep" << endl;
    }
    else if (linklist.qmf==2) {
        out.eprofile(linklist);
        cout << "printed first timestep" << endl;
        exit(1);
    }
    else if (linklist.qmf==4) {
        out.eccout(linklist);
        cout << "eccentricity printed" << endl;
    }


//    out.eprofile(linklist);
    while ((linklist.t<linklist.tend)&&(linklist.number_part<linklist.n())) {

        linklist.cfon=1;
        // adams<2>(dt,cc,&idealhydro3<2>,linklist);  // quicker method, not as accurate
        rungeKutta2<2>(dt,&idealhydro3<2>,linklist);  // slower method, more accurate

        if (linklist.cf>0) out.FOprint(linklist); // prints out frozen out SPH particles





        //if you add more points to print off you must also change LinkList<D>::setup and multiply steps=floor(tend-t0)+1; by the extra number of print offs / 1fm/c

        if (linklist.qmf==3) {
            double tsub=linklist.t-floor(linklist.t);
            if (tsub<(0.0+dt*0.99)||(tsub>=1-+dt*0.99)) // uncomment if you want to observe energydensity profile, conservation of energy or do a Gubser check
            {
                cout << "t=" << linklist.t <<endl;  // outputs time step
                out.eprofile(linklist);   // energy density profile
//        out.conservation(linklist); // conservation of energy
//        out.gubcheckux(linklist); // gubser test
//        out.gubcheckuy(linklist); // gubser test

            }
            else if ((tsub<(0.5+dt*0.5))&&(tsub>=(0.5-+dt*0.5))) // uncomment if you want to observe energydensity profile, conservation of energy or do a Gubser check
            {
                cout << "t=" <<  linklist.t <<endl;  // outputs time step
                out.eprofile(linklist);   // energy density profile
//        out.conservation(linklist); // conservation of energy
//        out.gubcheckux(linklist); // gubser test
//        out.gubcheckuy(linklist); // gubser test

            }

        }


    }

//    out.conservation(linklist);


    cout << "done" << endl;
    linklist.endEV();
}

void vSimulation(double dt,LinkList<2> &linklist)
{
    cout << "Ready to start hydrodynamics\n";




//     int cc=0; // uncomment for adams method

    linklist.frzc=0;
    linklist.cf=0;

    Output<2> out(linklist);

    linklist.t=linklist.t0;

    if (linklist.qmf==1) {
        out.eprofile(linklist);
        cout << "printed first timestep" << endl;
    }
    else if(linklist.qmf==4) {
        out.eccout(linklist);
        cout << "eccentricity printed" << endl;
    }

    linklist.Ez=0;
    while ((linklist.t<linklist.tend)&&(linklist.number_part<linklist.n())) {
        linklist.cfon=1;

        //vadams<2>(dt,cc,&vischydro3<2>,linklist); // quicker method, not as accurate
        vrungeKutta2<2>(dt,&vischydro3<2>,linklist);  // slower method, more accurate

        if (linklist.cf>0)
        {
            out.vFOprint(linklist); // prints off list of frozen out SPH particles
            // out.SCprint(linklist); // uncoment to print off entropy conservation
        }
//         double tsub=linklist.t-floor(linklist.t);
        //if you add more points to print, must also change LinkList<D>::setup and multiply steps=floor(tend-t0)+1; by the extra number of print offs / 1fm/c
//         if (tsub<(0.0+dt*0.99)||(tsub>=1-+dt*0.99)) // uncomment if you want to observe energydensity profile, conservation of energy or do a Gubser check
//         {
//        cout << "t=" << linklist.t <<endl;  // outputs time step
//             out.eprofile(linklist);   // energy density profile
//        out.conservation(linklist); // conservation of energy

//         }


    }
//    out.conservation(linklist);


    cout << "done" << endl;
    linklist.endEV();

}

void svSimulation(double dt,LinkList<2> &linklist)
{
    cout << "Ready to start hydrodynamics\n";




//     int cc=0;    //unused variable
    linklist.frzc=0;
    linklist.cf=0;


    Output<2> out(linklist);

    BBMG<2> bbmg(linklist);
    bbmg.initial(linklist);
    cout << "started bbmg" << endl;

    linklist.t=linklist.t0;
//    out.sveprofile(linklist);

    if (linklist.qmf==1||linklist.qmf==3) {
        out.sveprofile(linklist);
        cout << "printed first timestep" << endl;
        linklist.conservation_entropy();
        cout << "t=" << linklist.t << " S=" << linklist.S << endl;
        if (linklist.qmf==1) exit(0);
    }
    else if(linklist.qmf==4) {
        out.eccout(linklist);
        cout << "eccentricity printed" << endl;
        exit(0);
    }

    //out.sveprofile(linklist);
    linklist.Ez=0;
    while ((linklist.t<linklist.tend)&&(linklist.number_part<linklist.n())) {
        linklist.cfon=1;


        svrungeKutta2<2>(dt,&shear<2>,linklist);
        linklist.conservation_entropy();
        cout << "t=" << linklist.t << " " <<  linklist.Eloss << " " << linklist.S <<  endl;
        out.sveprofile(linklist);
//        out.sveprofile(linklist);
//        bbmg.propogate(linklist);


        if (linklist.cf>0) out.svFOprint(linklist);

        if (linklist.qmf==3) {
            double tsub=linklist.t-floor(linklist.t);
            // if you add more points to print, must also change LinkList<D>::setup and multiply steps=floor(tend-t0)+1; by the extra number of print offs / 1fm/c
            if (tsub<(0.0+dt*0.99)||(tsub>=1-+dt*0.99)) // uncomment if you want to observe energydensity profile, conservation of energy or do a Gubser check
            {
                linklist.conservation_entropy();
                cout << "t=" << linklist.t << " S=" << linklist.S << endl;  // outputs time step
                out.sveprofile(linklist);   // energy density profile
                cout << "eloss= " << linklist.t << " " <<  linklist.Eloss << endl;
                //out.conservation(linklist); // conservation of energy

            }
//         else if ((tsub<(0.2+dt*0.5))&&(tsub>=(0.2-+dt*0.5))) // uncomment if you want to observe energydensity profile, conservation of energy or do a Gubser check
//         {
//        cout << "t=" <<  linklist.t <<endl;  // outputs time step
//             out.sveprofile(linklist);   // energy density profile
////        out.conservation(linklist); // conservation of energy
////        out.gubcheckux(linklist); // gubser test
////        out.gubcheckuy(linklist); // gubser test
//
//         }
//         else if ((tsub<(0.4+dt*0.5))&&(tsub>=(0.4-+dt*0.5))) // uncomment if you want to observe energydensity profile, conservation of energy or do a Gubser check
//         {
//        cout << "t=" <<  linklist.t <<endl;  // outputs time step
//             out.sveprofile(linklist);   // energy density profile
////        out.conservation(linklist); // conservation of energy
////        out.gubcheckux(linklist); // gubser test
////        out.gubcheckuy(linklist); // gubser test
//
//         }
            else if ((tsub<(0.5+dt*0.5))&&(tsub>=(0.5-+dt*0.5))) // uncomment if you want to observe energydensity profile, conservation of energy or do a Gubser check
            {
                cout << "t=" <<  linklist.t <<endl;  // outputs time step
                out.sveprofile(linklist);   // energy density profile
//        out.conservation(linklist); // conservation of energy
//        out.gubcheckux(linklist); // gubser test
//        out.gubcheckuy(linklist); // gubser test

            }
//         else if ((tsub<(0.8+dt*0.5))&&(tsub>=(0.8-+dt*0.5))) // uncomment if you want to observe energydensity profile, conservation of energy or do a Gubser check
//         {
//        cout << "t=" <<  linklist.t <<endl;  // outputs time step
//             out.sveprofile(linklist);   // energy density profile
////        out.conservation(linklist); // conservation of energy
////        out.gubcheckux(linklist); // gubser test
////        out.gubcheckuy(linklist); // gubser test
//
//         }

        }

    }
    //out.conservation(linklist,pcount);


    linklist.endEV();
}

void BSQSimulation(double dt,LinkList<2> &linklist)
{
    cout << "Ready to start hydrodynamics\n";




//     int cc=0;  ///unneded variable
    linklist.frzc=0;
    linklist.cf=0;


    Output<2> out(linklist);

    BBMG<2> bbmg(linklist);
    bbmg.initial(linklist);
    cout << "started bbmg" << endl;

    linklist.t=linklist.t0;
//    out.sveprofile(linklist);

    if (linklist.qmf==1||linklist.qmf==3) {
        out.bsqsveprofile(linklist);
        cout << "printed first timestep" << endl;
        linklist.conservation_entropy();
        cout << "t=" << linklist.t << " S=" << linklist.S << endl;
        if (linklist.qmf==1) exit(0);
    }
    else if(linklist.qmf==4) {
        out.eccout(linklist);
        cout << "eccentricity printed" << endl;
        exit(0);
    }


	cout << "Now let's do the main evolution!" << endl;
    //out.sveprofile(linklist);
    linklist.Ez=0;
    while ((linklist.t<linklist.tend)&&(linklist.number_part<linklist.n())) {
        linklist.cfon=1;


	cout << "Entering here:" << endl;
        bsqrungeKutta2<2>(dt,&BSQshear<2>,linklist);
        linklist.conservation_entropy();
        cout << "t=" << linklist.t << " " <<  linklist.Eloss << " " << linklist.S <<  endl;
        out.bsqsveprofile(linklist);
//        out.sveprofile(linklist);
//        bbmg.propogate(linklist);


        if (linklist.cf>0) out.bsqsvFOprint(linklist);

        if (linklist.qmf==3) {
            double tsub=linklist.t-floor(linklist.t);
            // if you add more points to print, must also change LinkList<D>::setup and multiply steps=floor(tend-t0)+1; by the extra number of print offs / 1fm/c
            if (tsub<(0.0+dt*0.99)||(tsub>=1-+dt*0.99)) // uncomment if you want to observe energydensity profile, conservation of energy or do a Gubser check
            {
                linklist.conservation_entropy();
                cout << "t=" << linklist.t << " S=" << linklist.S << endl;  // outputs time step
                out.bsqsveprofile(linklist);   // energy density profile
                cout << "eloss= " << linklist.t << " " <<  linklist.Eloss << endl;
                //out.conservation(linklist); // conservation of energy

            }
//         else if ((tsub<(0.2+dt*0.5))&&(tsub>=(0.2-+dt*0.5))) // uncomment if you want to observe energydensity profile, conservation of energy or do a Gubser check
//         {
//        cout << "t=" <<  linklist.t <<endl;  // outputs time step
//             out.sveprofile(linklist);   // energy density profile
////        out.conservation(linklist); // conservation of energy
////        out.gubcheckux(linklist); // gubser test
////        out.gubcheckuy(linklist); // gubser test
//
//         }
//         else if ((tsub<(0.4+dt*0.5))&&(tsub>=(0.4-+dt*0.5))) // uncomment if you want to observe energydensity profile, conservation of energy or do a Gubser check
//         {
//        cout << "t=" <<  linklist.t <<endl;  // outputs time step
//             out.sveprofile(linklist);   // energy density profile
////        out.conservation(linklist); // conservation of energy
////        out.gubcheckux(linklist); // gubser test
////        out.gubcheckuy(linklist); // gubser test
//
//         }
            else if ((tsub<(0.5+dt*0.5))&&(tsub>=(0.5-+dt*0.5))) // uncomment if you want to observe energydensity profile, conservation of energy or do a Gubser check
            {
                cout << "t=" <<  linklist.t <<endl;  // outputs time step
                out.bsqsveprofile(linklist);   // energy density profile
//        out.conservation(linklist); // conservation of energy
//        out.gubcheckux(linklist); // gubser test
//        out.gubcheckuy(linklist); // gubser test

            }
//         else if ((tsub<(0.8+dt*0.5))&&(tsub>=(0.8-+dt*0.5))) // uncomment if you want to observe energydensity profile, conservation of energy or do a Gubser check
//         {
//        cout << "t=" <<  linklist.t <<endl;  // outputs time step
//             out.sveprofile(linklist);   // energy density profile
////        out.conservation(linklist); // conservation of energy
////        out.gubcheckux(linklist); // gubser test
////        out.gubcheckuy(linklist); // gubser test
//
//         }

        }

    }
    //out.conservation(linklist,pcount);


    linklist.endEV();
}

void Simulation(double dt,LinkList<3> &linklist)  // setup for 3+1, equation of motion not complete!!!
{
    double t=0.6;
    int cc=0;


    Output<3> out(linklist);


    linklist.t=t;
    idealhydro3<3>(linklist) ;
    out.eprofile(linklist);



    linklist.Ez=0;
    while (t<linklist.tend) {
        adams<3>(dt,cc,&idealhydro3<3>,linklist);



        double tsub=linklist.t-floor(linklist.t);
        if ((tsub>0.599)&&(tsub<0.609))
        {

            out.eprofile(linklist);
            out.conservation(linklist);
            out.gubcheckux(linklist);
            out.gubcheckuy(linklist);

        }

    }
    linklist.endEV();

}

template <int D>
void idealhydro3(LinkList<D>  &linklist) // ideal Equations of motion, only set up completely for 2+1 at the moment
{


    linklist.initiate();  // initiates linklist


    for(int i=0; i<linklist.n(); i++)
    {
        linklist.optimization(i); // calculates entropy density
    }

    int curfrz=0;
    for(int i=0; i<linklist.n(); i++)
    {
        //  Computes gamma and velocity

        linklist._p[i].calc(linklist.t); // calculates gamma, flow vectors, and updates EOS

        linklist._p[i].returnA(); // returns A form EOM
//            if(linklist._p[i].EOSs() < 0)
        if (linklist.cfon==1) linklist._p[i].frzcheck(linklist.t,curfrz,linklist.n()); // checks if SPH particles have frozen out


    }

    if (linklist.cfon==1) // creates list of frozen out SPH particles
    {
        linklist.number_part+=curfrz;
        linklist.list.resize(curfrz);
    }
    if (linklist.rk2==1) linklist.conservation();    // calculates conservation of energy, only on first Runge kutta step

    int m=0;
    for(int i=0; i<linklist.n(); i++)
    {
        //      Computes gradients to obtain dsigma/dt
        linklist.optimization2(i);
        linklist._p[i].sigset(linklist.t); // returns another value for EOM
        if ((linklist._p[i].Freeze==3)&&(linklist.cfon==1)) // if SPH particles have remained frozen out after two time steps, allows to freezeout
        {
            linklist.list[m]=i;
            linklist._p[i].Freeze=4;
            ++m;
        }

    }

    linklist.conservation_Ez(); // calculates dEz for conservation of Energy
    if (linklist.first==1) //checks if t=t0, if so sets the initial energy E0
    {
        linklist.conservation_entropy();
        linklist.conservation_E();
    }



    //calculate matrix elements
    for(int i=0; i<linklist.n(); i++)
    {

        // set the Mass and the Force matrix
        double *M,*F;
        M=new double[D*D];
        F=new double[D];
        M[0]=linklist._p[i].Agam*linklist._p[i].u.x[0]*linklist._p[i].u.x[0]+linklist._p[i].EOSw()*linklist._p[i].gamma;
        M[3]=linklist._p[i].Agam*linklist._p[i].u.x[1]*linklist._p[i].u.x[1]+linklist._p[i].EOSw()*linklist._p[i].gamma;
        M[1]=linklist._p[i].Agam*linklist._p[i].u.x[0]*linklist._p[i].u.x[1];

        F[0]=linklist._p[i].Agam2*linklist._p[i].u.x[0]-linklist._p[i].gradP.x[0];
        F[1]=linklist._p[i].Agam2*linklist._p[i].u.x[1]-linklist._p[i].gradP.x[1];




        double det=M[0]*M[3]-M[1]*M[1]; // inverts matrix
        double MI[4];
        MI[0]=M[3]/det;
        MI[1]=-M[1]/det;
        MI[3]=M[0]/det;
        linklist._p[i].du_dt.x[0]=F[0]*MI[0]+F[1]*MI[1];
        linklist._p[i].du_dt.x[1]=F[0]*MI[1]+F[1]*MI[3];



        linklist._p[i].div_u = (1./ linklist._p[i].gamma)*inner( linklist._p[i].u, linklist._p[i].du_dt) - ( linklist._p[i].gamma/ linklist._p[i].sigma)* linklist._p[i].dsigma_dt ;

        delete [] M;
        delete [] F;

    }
    if (linklist.cfon==1) linklist.freezeout(curfrz); // calculates normals from freezeout hypersurface





    linklist.destroy();
}

template <int D>
void vischydro3(LinkList<D>  &linklist)  // bulk Equations of motion, only set up completely for 2+1 at the moment
{
    linklist.initiate();
    for(int i=0; i<linklist.n(); i++)
    {
        linklist.voptimization(i);
        if ((linklist._p[i].eta<0)||isnan(linklist._p[i].eta))
        {
            cout << "neg entropy " <<  linklist._p[i].EOST()*197.3   << " " << linklist._p[i].eta << endl;


            //double seta = 0;
            Vector<int,D> jj;
            int b;
            for(jj.x[0]=-2; jj.x[0]<=2; jj.x[0]++)
            {
                for(jj.x[1]=-2; jj.x[1]<=2; jj.x[1]++)
                {

                    b=linklist.lead[linklist.triToSum(linklist.dael[i]+jj, linklist.size)];
                    while(b!=-1 )
                    {
                        if (linklist._p[b].eta_sigma<0)
                        {
                            cout << b << " "   <<linklist._p[b].eta_sigma<< " "  <<  linklist._p[i].EOST()*197.3  << " " << linklist._p[i].Bulk << endl;
                            cout << linklist._p[b].tauRelax << " "  <<  linklist._p[i].zeta << endl;
                        }
                        b=linklist.link[b];
                    }
                }
            }

            exit(1);
        }

    }
    int curfrz=0;
    for(int i=0; i<linklist.n(); i++)
    {
        //  Computes gamma and velocity

        linklist._p[i].calc(linklist.t);
        linklist._p[i].setvisc(linklist.etaconst,linklist.bvf,linklist.svf,linklist.zTc,linklist.sTc,linklist.zwidth,linklist.visc);
        if (linklist.cfon==1) linklist._p[i].frzcheck(linklist.t,curfrz,linklist.n());
    }
    if (linklist.cfon==1)
    {
        linklist.number_part+=curfrz;
        linklist.list.resize(curfrz);
    }
    int m=0;
    for(int i=0; i<linklist.n(); i++)
    {
        //      Computes gradients to obtain dsigma/dt
        linklist.voptimization2(i,linklist.t);
        linklist._p[i].vsigset(linklist.t);
        if ((linklist._p[i].Freeze==3)&&(linklist.cfon==1))
        {
            linklist.list[m]=i;
            linklist._p[i].Freeze=4;
            ++m;
        }

    }

    if (linklist.rk2==1) linklist.vconservation();
    linklist.vconservation_Ez();



    //calculate matrix elements
    for(int i=0; i<linklist.n(); i++)
    {


        // set the Mass and the Force
        double *M,*F;
        M=new double[D*D];
        F=new double[D];
        M[0]=linklist._p[i].Agam*linklist._p[i].u.x[0]*linklist._p[i].u.x[0]+linklist._p[i].C*linklist._p[i].gamma;
        M[3]=linklist._p[i].Agam*linklist._p[i].u.x[1]*linklist._p[i].u.x[1]+linklist._p[i].C*linklist._p[i].gamma;
        M[1]=linklist._p[i].Agam*linklist._p[i].u.x[0]*linklist._p[i].u.x[1];

        F[0]=linklist._p[i].Agam2*linklist._p[i].u.x[0]-(linklist._p[i].gradP.x[0]+linklist._p[i].gradBulk.x[0]);
        F[1]=linklist._p[i].Agam2*linklist._p[i].u.x[1]-(linklist._p[i].gradP.x[1]+linklist._p[i].gradBulk.x[1]);



        double MI[4];

        double det=M[0]*M[3]-M[1]*M[1];
        MI[0]=M[3]/det;
        MI[1]=-M[1]/det;
        MI[3]=M[0]/det;
        linklist._p[i].du_dt.x[0]=F[0]*MI[0]+F[1]*MI[1];
        linklist._p[i].du_dt.x[1]=F[0]*MI[1]+F[1]*MI[3];

        linklist._p[i].div_u = (1./ linklist._p[i].gamma)*inner( linklist._p[i].u, linklist._p[i].du_dt) - ( linklist._p[i].gamma/ linklist._p[i].sigma)* linklist._p[i].dsigma_dt ;

        linklist._p[i].bigtheta=linklist._p[i].div_u*linklist.t+linklist._p[i].gamma;

        linklist._p[i].detasigma_dt = -linklist._p[i].bigPI/linklist._p[i].EOST()*linklist._p[i].bigtheta/linklist._p[i].sigma;

        linklist._p[i].dBulk_dt = (-linklist._p[i].zeta/linklist._p[i].sigma*linklist._p[i].bigtheta - linklist._p[i].Bulk/linklist._p[i].gamma )/linklist._p[i].tauRelax;

        delete [] M;
        delete [] F;

    }

    if (linklist.cfon==1) linklist.vfreezeout(curfrz);



    linklist.destroy();
}

template <int D>
void shear(LinkList<D>  &linklist)  // shear+bulk Equations of motion, only set up completely for 2+1 at the moment
{



    linklist.setshear();
    linklist.initiate();




    for(int i=0; i<linklist.n(); i++)
    {
        linklist.optimization(i);


        if ((linklist._p[i].eta<0)||isnan(linklist._p[i].eta))
        {
            cout << i <<  " neg entropy " <<  linklist._p[i].EOST()*197.3   << " " << linklist._p[i].eta << endl;

            linklist._p[i].eta=0;

//             int a,b;
//             a=i;
//            double sigma = 0.;
//                Vector<int,D> ii;
//            for(ii.x[0]=-2; ii.x[0]<=2; ii.x[0]++)
//                {
//                for(ii.x[1]=-2; ii.x[1]<=2; ii.x[1]++)
//                {

//                     b=linklist.lead[linklist.triToSum(linklist.dael[a]+ii, linklist.size)];
//                     while(b!=-1 )
//                     {
//                     double kern=linklist.kernel(linklist._p[b].r-linklist._p[a].r);
//                        sigma = sigma + linklist._p[b].sigmaweight*kern;
//
//
//                       b=linklist.link[b];
//                       }
//                }
//            }
//            cout << "final entropy=" << sigma << " eta_sigma=" << linklist._p[a].eta_sigma << " r=" << linklist._p[a].r << endl;
//
//            exit(1);
        }

    }


    int curfrz=0;
    for(int i=0; i<linklist.n(); i++)
    {
        //  Computes gamma and velocity

        linklist._p[i].calc(linklist.t);
        linklist._p[i].setvisc(linklist.etaconst,linklist.bvf,linklist.svf,linklist.zTc,linklist.sTc,linklist.zwidth,linklist.visc);
        if (linklist.cfon==1) linklist._p[i].frzcheck(linklist.t,curfrz,linklist.n());

    }
    if (linklist.cfon==1)
    {
        linklist.number_part+=curfrz;
        linklist.list.resize(curfrz);
    }

    int m=0;
    for(int i=0; i<linklist.n(); i++)
    {
        //      Computes gradients to obtain dsigma/dt
        linklist.svoptimization2(i,linklist.t,curfrz);
        linklist._p[i].dsigma_dt = -linklist._p[i].sigma*(linklist._p[i].gradV.x[0][0]+linklist._p[i].gradV.x[1][1]) ;

        linklist._p[i].svsigset(linklist.t,i);
        if ((linklist._p[i].Freeze==3)&&(linklist.cfon==1))
        {
            linklist.list[m]=i;
            linklist._p[i].Freeze=4;
            ++m;
        }

    }

    if (linklist.rk2==1) linklist.svconservation();
    linklist.svconservation_Ez();


    //calculate matrix elements
    for(int i=0; i<linklist.n(); i++)
    {
        double gamt=1./linklist._p[i].gamma/linklist._p[i].stauRelax;
        double pre=linklist._p[i].eta_o_tau/2./linklist._p[i].gamma;
        //p4=gamt-linklist._p[i].sigl*4./3.;
        double p1=gamt-4./3./linklist._p[i].sigma*linklist._p[i].dsigma_dt+1./linklist.t/3.;
        Vector<double,D>  minshv=rowp1(0,linklist._p[i].shv);
        //p2=linklist._p[i].setas*gamt;
        Matrix <double,D,D> partU=linklist._p[i].gradU+transpose(linklist._p[i].gradU);

        // set the Mass and the Force
        Matrix <double,D,D> M=linklist._p[i].Msub(i);
        Vector<double,D> F=linklist._p[i].Btot*linklist._p[i].u+ linklist._p[i].gradshear -(linklist._p[i].gradP+linklist._p[i].gradBulk+ linklist._p[i].divshear);
        // shear contribution
        F+=pre*linklist._p[i].v*partU+p1*minshv;


        double det=deter(M);
        Matrix <double,D,D> MI;
        MI.x[0][0]=M.x[1][1]/det;
        MI.x[0][1]=-M.x[0][1]/det;
        MI.x[1][0]=-M.x[1][0]/det;
        MI.x[1][1]=M.x[0][0]/det;
        linklist._p[i].du_dt.x[0]=(F.x[0]*MI.x[0][0]+F.x[1]*MI.x[0][1]);
        linklist._p[i].du_dt.x[1]=(F.x[0]*MI.x[1][0]+F.x[1]*MI.x[1][1]);




        Matrix <double,D,D> ulpi=linklist._p[i].u*colp1(0,linklist._p[i].shv);

        double vduk=inner(linklist._p[i].v,linklist._p[i].du_dt);

        Matrix <double,D,D> Ipi=-linklist._p[i].eta_o_tau/3.*(linklist._p[i].Imat+linklist._p[i].uu)+4./3.*linklist._p[i].pimin;

        linklist._p[i].div_u = (1./ linklist._p[i].gamma)*inner( linklist._p[i].u, linklist._p[i].du_dt) - ( linklist._p[i].gamma/ linklist._p[i].sigma)* linklist._p[i].dsigma_dt ;
        linklist._p[i].bigtheta=linklist._p[i].div_u*linklist.t+linklist._p[i].gamma;


        Matrix <double,D,D> sub=linklist._p[i].pimin+linklist._p[i].shv.x[0][0]/linklist._p[i].g2*linklist._p[i].uu-1./linklist._p[i].gamma*linklist._p[i].piutot;


        linklist._p[i].inside=linklist.t*(inner((-minshv+linklist._p[i].shv.x[0][0]*linklist._p[i].v),linklist._p[i].du_dt)- con2(sub,linklist._p[i].gradU)    -linklist._p[i].gamma*linklist.t*linklist._p[i].shv33);
        linklist._p[i].detasigma_dt =1./linklist._p[i].sigma/linklist._p[i].EOST()*( -linklist._p[i].bigPI*linklist._p[i].bigtheta+linklist._p[i].inside);




        linklist._p[i].dBulk_dt = (-linklist._p[i].zeta/linklist._p[i].sigma*linklist._p[i].bigtheta - linklist._p[i].Bulk/linklist._p[i].gamma )/linklist._p[i].tauRelax;

        Matrix <double,D,D> ududt=linklist._p[i].u*linklist._p[i].du_dt;


        linklist._p[i].dshv_dt= -gamt*(linklist._p[i].pimin+linklist._p[i].setas*0.5*partU)-0.5*linklist._p[i].eta_o_tau*(ududt+transpose(ududt))+linklist._p[i].dpidtsub()-vduk*(ulpi+transpose(ulpi)+(1/linklist._p[i].gamma)*Ipi)+linklist._p[i].sigl*Ipi;



    }


    if (linklist.cfon==1) linklist.svfreezeout(curfrz);


    linklist.destroy();
}

template <int D>
void BSQshear(LinkList<D>  &linklist)  // shear+bulk Equations of motion, only set up completely for 2+1 at the moment
{



    linklist.setshear();
    linklist.initiate();




    for(int i=0; i<linklist.n(); i++)
    {
		cout << "Entering this loop: i = " << i << endl;
        int curfrz=0;//added by Christopher Plumberg to get compilation
        linklist.bsqsvoptimization(i);    // NOT bsqsvoptimization2!!! fix arguments accordingly!!!

        if ((linklist._p[i].eta<0)||isnan(linklist._p[i].eta))
        {
            cout << i <<  " neg entropy " <<  linklist._p[i].EOST()*197.3   << " " << linklist._p[i].eta << endl;

            linklist._p[i].eta=0;

//             int a,b;
//             a=i;
//            double sigma = 0.;
//                Vector<int,D> ii;
//            for(ii.x[0]=-2; ii.x[0]<=2; ii.x[0]++)
//                {
//                for(ii.x[1]=-2; ii.x[1]<=2; ii.x[1]++)
//                {

//                     b=linklist.lead[linklist.triToSum(linklist.dael[a]+ii, linklist.size)];
//                     while(b!=-1 )
//                     {
//                     double kern=linklist.kernel(linklist._p[b].r-linklist._p[a].r);
//                        sigma = sigma + linklist._p[b].sigmaweight*kern;
//
//
//                       b=linklist.link[b];
//                       }
//                }
//            }
//            cout << "final entropy=" << sigma << " eta_sigma=" << linklist._p[a].eta_sigma << " r=" << linklist._p[a].r << endl;
//
//            exit(1);
        }

    }


    int curfrz=0;
    for(int i=0; i<linklist.n(); i++)
    {
        //  Computes gamma and velocity

        linklist._p[i].calcbsq(linklist.t); //resets EOS!!
        linklist._p[i].setvisc(linklist.etaconst,linklist.bvf,linklist.svf,linklist.zTc,linklist.sTc,linklist.zwidth,linklist.visc);
        if (linklist.cfon==1) linklist._p[i].frzcheck(linklist.t,curfrz,linklist.n());

    }
    if (linklist.cfon==1)
    {
        linklist.number_part+=curfrz;
        linklist.list.resize(curfrz);
    }

    int m=0;
    for(int i=0; i<linklist.n(); i++)
    {
        //      Computes gradients to obtain dsigma/dt
        linklist.bsqsvoptimization2(i,linklist.t,curfrz);

        linklist._p[i].dsigma_dt = -linklist._p[i].sigma
									*( linklist._p[i].gradV.x[0][0]
										+ linklist._p[i].gradV.x[1][1]) ;
cout << "CHECK dsigma_dt: " << i << "   " << linklist._p[i].dsigma_dt
		<< "   " << linklist._p[i].sigma << "   " << linklist._p[i].gradV.x[0][0]
		<< "   " << linklist._p[i].gradV.x[1][1] << endl;

        linklist._p[i].bsqsvsigset(linklist.t,i);
        if ((linklist._p[i].Freeze==3)&&(linklist.cfon==1))
        {
            linklist.list[m]=i;
            linklist._p[i].Freeze=4;
            ++m;
        }

    }

    if (linklist.rk2==1) linklist.bsqsvconservation();
    linklist.bsqsvconservation_Ez();


    //calculate matrix elements
    for(int i=0; i<linklist.n(); i++)
    {
        double gamt=1./linklist._p[i].gamma/linklist._p[i].stauRelax;
        double pre=linklist._p[i].eta_o_tau/2./linklist._p[i].gamma;
        //p4=gamt-linklist._p[i].sigl*4./3.;
        double p1=gamt-4./3./linklist._p[i].sigma*linklist._p[i].dsigma_dt+1./linklist.t/3.;
        Vector<double,D>  minshv=rowp1(0,linklist._p[i].shv);
        //p2=linklist._p[i].setas*gamt;
        Matrix <double,D,D> partU=linklist._p[i].gradU+transpose(linklist._p[i].gradU);

        // set the Mass and the Force
        Matrix <double,D,D> M=linklist._p[i].Msub(i);
        Vector<double,D> F=linklist._p[i].Btot*linklist._p[i].u+ linklist._p[i].gradshear -(linklist._p[i].gradP+linklist._p[i].gradBulk+ linklist._p[i].divshear);
        // shear contribution
        F+=pre*linklist._p[i].v*partU+p1*minshv;


        double det=deter(M);
        Matrix <double,D,D> MI;
        MI.x[0][0]=M.x[1][1]/det;
        MI.x[0][1]=-M.x[0][1]/det;
        MI.x[1][0]=-M.x[1][0]/det;
        MI.x[1][1]=M.x[0][0]/det;
        linklist._p[i].du_dt.x[0]=F.x[0]*MI.x[0][0]+F.x[1]*MI.x[0][1];
        linklist._p[i].du_dt.x[1]=F.x[0]*MI.x[1][0]+F.x[1]*MI.x[1][1];




        Matrix <double,D,D> ulpi=linklist._p[i].u*colp1(0,linklist._p[i].shv);

        double vduk=inner(linklist._p[i].v,linklist._p[i].du_dt);

        Matrix <double,D,D> Ipi=-linklist._p[i].eta_o_tau/3.*(linklist._p[i].Imat+linklist._p[i].uu)+4./3.*linklist._p[i].pimin;

        linklist._p[i].div_u = (1./ linklist._p[i].gamma)*inner( linklist._p[i].u, linklist._p[i].du_dt)
								- ( linklist._p[i].gamma/ linklist._p[i].sigma)
									* linklist._p[i].dsigma_dt ;
        linklist._p[i].bigtheta=linklist._p[i].div_u*linklist.t+linklist._p[i].gamma;

cout << "CHECK div_u: " << i << "   " << linklist._p[i].div_u
		<< "   " << linklist._p[i].gamma
		<< "   " << linklist._p[i].u
		<< "   " << linklist._p[i].du_dt
		<< "   " << inner( linklist._p[i].u, linklist._p[i].du_dt)
		<< "   " << linklist._p[i].sigma << endl;
cout << "CHECK bigtheta: " << i << "   " << linklist._p[i].bigtheta
		<< "   " << linklist.t
		<< "   " << linklist._p[i].gamma << endl;

        Matrix <double,D,D> sub=linklist._p[i].pimin+linklist._p[i].shv.x[0][0]/linklist._p[i].g2*linklist._p[i].uu-1./linklist._p[i].gamma*linklist._p[i].piutot;


        linklist._p[i].inside=linklist.t*(inner((-minshv+linklist._p[i].shv.x[0][0]*linklist._p[i].v),linklist._p[i].du_dt)- con2(sub,linklist._p[i].gradU)    -      linklist._p[i].gamma*linklist.t*linklist._p[i].shv33);
        linklist._p[i].detasigma_dt =1./linklist._p[i].sigma/linklist._p[i].EOST()
										*( -linklist._p[i].bigPI*linklist._p[i].bigtheta
											+linklist._p[i].inside);
std::cout << "Check detasigma_dt: " << i << "   "
			<< linklist._p[i].sigma << "   "
			<< linklist._p[i].EOST() << "   "
			<< linklist._p[i].bigPI << "   "
			<< linklist._p[i].bigtheta << "   "
			<< linklist._p[i].inside << std::endl;




        linklist._p[i].dBulk_dt = (-linklist._p[i].zeta/linklist._p[i].sigma*linklist._p[i].bigtheta - linklist._p[i].Bulk/linklist._p[i].gamma )/linklist._p[i].tauRelax;

        Matrix <double,D,D> ududt=linklist._p[i].u*linklist._p[i].du_dt;


        linklist._p[i].dshv_dt= -gamt*(linklist._p[i].pimin+linklist._p[i].setas*0.5*partU)-0.5*linklist._p[i].eta_o_tau*(ududt+transpose(ududt))+linklist._p[i].dpidtsub()-vduk*(ulpi+transpose(ulpi)+(1/linklist._p[i].gamma)*Ipi)+linklist._p[i].sigl*Ipi;

        //linklist._p[i].drhoB_dt=-linklist._p[i].rhoB*linklist._p[i].sigma*linklist._p[i].bigtheta;
        //linklist._p[i].drhoS_dt=-linklist._p[i].rhoS*linklist._p[i].sigma*linklist._p[i].bigtheta;
        //linklist._p[i].drhoQ_dt=-linklist._p[i].rhoQ*linklist._p[i].sigma*linklist._p[i].bigtheta;

    }


    if (linklist.cfon==1) linklist.bsqsvfreezeout(curfrz);


    linklist.destroy();
}
