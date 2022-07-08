#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <string.h>

using namespace std;


#include "vector.h"
#include "tables.h"
#include "particle.h"
#include "enteric.h"
#include "hydrosim.h"
#include "LinkList.h"
#include "eostables.h"
#include "matrix.h"



char ifolder []="inputfiles/";
string ofolder ("./outputfiles/");
double freezeoutT;
double freezeoutB;
double freezeoutS;
double freezeoutQ;
double zconst=1/(4*PI);
eostableshigh *ETH;
eostableslow *ETL;
int nETH,nETL;
Matrix <double,2,2> Imat;
int hcor;

template <int D>
eos Particle<D>::EOS = eos();

// This is the main file of v-USPhydro.  The files are organized as follows:
// EOS- eos.cpp/eos.h but tables to read them in are in tables.h and files to read them in happen within enetric.h/enetric.cpp
// Hydro evolution equations- primarily hydrosim.cpp/hydrosim.h but parts are found withint subroutines of Linklist.h and particle.h
// Initial conditions- enetric.h/enetric.cpp but some parts are processed within Linklist.h
// Math framework for vectors and matrices: vector.h/vector.cpp and matrix.h/matrix.cpp
// Integration subroutines- rungekutta.h
// Linklist- Linklist.h
// Calculation of freezeout hypersurface - particles are indentified in hydrosim.cppp but normal vectors are calculated in Linklist.h
// output- output.h


// v-USPhydro is an SPH code that solves 2+1 relativistic ideal or viscous hydrodynamical equations.  In this version there is infrastructure to write the 3+1 version but it is not yet set up.  The Cooper Frye freeze-out and decays are found either in subfolders df (analytical solutions, no decays) or sampling (Monte Carlo sampling+decays).  The respective flow harmonics are calculated within each routine.

// Last updated JNH 19.8.2014

int main (int argc, char *argv[])
{
// initialize variables, reads in positions and velocities

    _inputIC ics;

    if (argv[1])
    {
        ics.man= argv[1];

        if (argv[2])
        {
            ics.rnum=argv[2];



            if (argv[3]&&argv[4])
            {
                stringstream s,s1;
                s << argv[3];
                s >> ics.start;

                s1 << argv[4];
                s1 >> ics.end;

//    cout << "Start at event " << ics.start  << " end at " << ics.end << endl;
                ics.on=1;
            }
            else if (argv[3])
            {
                stringstream s,s1;
                s << argv[3];
                s >> ics.start;

                ics.end=ics.start;

//    cout << "Run event " << ics.start  << endl;
                ics.on=1;
            }
            else {
                ics.on=0;
            }
        }

    }
    else
    {
        ics.man="manual2.dat";
        ics.on=0;
    }


    int di;
// double dt;
    Imat.identity();





//  cout << "2+1 v-USPhydro in C++\n";
//  cout << "by JNH\n";
    cout << "*******************************\n";

//      Include file that has h,timestep,dimensions and files for EOS/IC's
//      File should be stores in the inputfile directory

    di=2; // define number of dimensions HERE

    if (di==2) // two dimensions
    {
        LinkList<2> linklist;
        manualenter<2>(ics,linklist);
        if (linklist.visc==0) // ideal
        {
            Simulation(ics.dt,linklist);
            if (linklist.fcount!=0)
            {

                for (int i=1; i<linklist.fcount; i++) // runs over other events
                {

                    linklist.cevent=i;
                    nextevent<2>(i, linklist);
                    cout << "Event " << linklist.fnum << endl;
                    Simulation(ics.dt,linklist);
                }

                delete [] linklist.filenames;
            }
        }
        else if (linklist.visc==1) // bulk
        {

            vSimulation(ics.dt,linklist);
            if (linklist.fcount!=0)
            {

                for (int i=1; i<linklist.fcount; i++) // runs over other events
                {

                    linklist.cevent=i;
                    nextevent<2>(i, linklist);
                    cout << "Event " << linklist.fnum << endl;
                    vSimulation(ics.dt,linklist);
                }
                delete [] linklist.filenames;

            }
        }
        else if (linklist.visc==3)  // shear+bulk
        {
            svSimulation(ics.dt,linklist);
            if (linklist.fcount!=0)
            {

                for (int i=1; i<linklist.fcount; i++) // runs over other events
                {

                    linklist.cevent=i;
                    nextevent<2>(i, linklist);
                    cout << "Event " << linklist.fnum << endl;
                    svSimulation(ics.dt,linklist);
                }

                delete [] linklist.filenames;
            }
        }
        else if (linklist.visc==4)  // shear+bulk+BSQ
        {
            BSQSimulation(ics.dt,linklist);
            if (linklist.fcount!=0)
            {

                for (int i=1; i<linklist.fcount; i++) // runs over other events
                {

                    linklist.cevent=i;
                    nextevent<2>(i, linklist);
                    cout << "Event " << linklist.fnum << endl;
                    BSQSimulation(ics.dt,linklist);
                }

                delete [] linklist.filenames;
            }
        }
        else
        {
            cout << "Error: visc not set up!" << endl;
        }
        string table ("table");
        if (linklist.eost==table)
        {
            delete [] ETH;
            delete [] ETL;
        }

    }
    else if (di==3) // three dimensions
    {
        LinkList<3> linklist;
        manualenter<3>(ics,linklist);
        Simulation(ics.dt,linklist);

        if (linklist.fcount!=0)
        {

            for (int i=1; i<linklist.fcount; i++)
            {

                linklist.cevent=i;
                nextevent<3>(i, linklist);
                cout << "Event " << linklist.fnum << endl;
                Simulation(ics.dt,linklist);
            }
            delete [] linklist.filenames;

        }

        string table ("table");
        if (linklist.eost==table)
        {
            delete [] ETH;
            delete [] ETL;
        }
    }



    return 0;
}
