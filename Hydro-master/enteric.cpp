#ifndef _ENTERIC_CPP_
#define _ENTERIC_CPP_

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <vector>
#include "mathdef.h"
#include "vector.h"
#include "tables.h"
#include "particle.h"
#include "enteric.h"
#include<dirent.h>
#include <stdlib.h>
#include <algorithm>
#include <cmath>

using namespace std;


void readEOStable()
{
    /////////////////////////////////////////////////////
    ////       Start reading EoS table           /////
    /////////////////////////////////////////////////////

    cout << "The EOS should be contained in two files (either .txt or .dat files).  The first contains the Temperature (in MeV's) and the first line should include the number of steps.\n";
    cout << "Enter file name of the first file:\n";
    cin.clear();
    cin.ignore (cin.rdbuf( )->in_avail( ));
    string filen;
    getline (cin,filen);

    readEOS_T(filen);

    cout << "The second EOS (either .txt or .dat files) contains the  Energy Density, Pressure, Entropy,   and c_s^2 (in MeV's) and the first line should include the number of steps.\n";
    cout << "Enter file name of the first file:\n";
    cin.clear();
    cin.ignore (cin.rdbuf( )->in_avail( ));
    string filen2;
    getline (cin,filen2);

    readEOS_p(filen2);

    /////////////////////////////////////////////////////
    ////       Finished reading EoS table           /////
    /////////////////////////////////////////////////////
}


void readEOS_T(string &firstry)
{
    string filename;
    filename = ifolder+firstry;
    //nameenter:

    FILE * myfile = fopen (filename.c_str(),"r");

    if ((!firstry.empty())&&(myfile!= NULL)) {
        fscanf(myfile,"%i \n",&nETH);

        ETH=new eostableshigh[nETH];
        int i=0;
        while (i<nETH ) {
            fscanf(myfile,"%lf     %*f  %*f \n",&ETH[i].T); //assumes T is in the formate MeV

            ETH[i].T/=197.3;
            ++i;
        }
        fclose(myfile);
        // cout << filename.c_str() << ": Input sucessful!\n";
    }
    else {
        cout << "Error: " << firstry.c_str() << " does not exist.  Please enter new file name\n";

        exit(1);
    }


}


void readEOS_p(string &firstry)
{

    string filename;
    filename = ifolder+firstry;
    //nameenter:
    FILE * myfile = fopen (filename.c_str(),"r");
    if ((!firstry.empty())&&(myfile!= NULL)) {
        fscanf(myfile,"%*i \n");
        int i=0;
        while (i<nETH ) {
            fscanf(myfile,"%lf     %lf    %lf    %lf \n",&ETH[i].e,&ETH[i].p,&ETH[i].s,&ETH[i].dtds); //assumes that e,p is in the format GeV/fm^3 and s in 1/fm^3
            ETH[i].e/=0.1973;
            ETH[i].p/=0.1973;
            ++i;
        }
        fclose(myfile);
        // cout << filename.c_str() << ": Input sucessful!\n";
    }
    else {
        cout << "Error: " << firstry.c_str() << " does not exist.  Please enter new file name\n";

        exit(1);
    }
}

void readEOS_lowS(string &firstry)
{


    string filename;
    filename = ifolder+firstry;
    //nameenter:
    FILE * myfile = fopen (filename.c_str(),"r");

    if ((!firstry.empty())&&(myfile!= NULL)) {
        fscanf(myfile,"%i \n",&nETL);
        ETL = new eostableslow[nETL];
        int i=0;
        while (i<nETL ) {
            fscanf(myfile,"%lf     %lf  %lf    %lf   %lf\n",&ETL[i].T,&ETL[i].e,&ETL[i].p,&ETL[i].s,&ETL[i].dtds); //assumes T is in the formate MeV
            ETL[i].T/=197.3;
            i++;
        }
        fclose(myfile);
        //  cout << filename.c_str() << ": Input sucessful!\n";
    }
    else {
        cout << "Error: " << firstry.c_str() << " does not exist.  Please enter new file name\n";

        exit(1);
    }


}


//2 dimensions

// rone's version
void readICs_gen(string &firstry, string &firstry2, int &_Ntable3,Particle<2> *&_p,double const& efcheck, int & numpart)
{
    string filename;

    filename = ifolder+firstry;
    //nameenter:
    FILE * myfile = fopen (filename.c_str(),"r");
    if ((!firstry.empty())&&(myfile!= NULL)) {
        fscanf(myfile,"%i \n",&_Ntable3);
        _p= new Particle<2>[_Ntable3];
        int i=0;

        while (i<=_Ntable3 ) {
            fscanf(myfile,"%lf    %lf    %lf\n",&_p[i].r.x[0],&_p[i].r.x[1],&_p[i].sigmaweight);
            _p[i].u.x[0]=0;
            _p[i].u.x[1]=0;
            _p[i].eta_sigma  = 1.;
            _p[i].Bulk = 0.;

            if (_p[i].e_sub>efcheck)
            {
                _p[i].Freeze=0;
            }
            else
            {
                _p[i].Freeze=4;
                ++numpart;
            }

            ++i;
        }


        fclose(myfile);
        cout << filename.c_str() << ": Input sucessful!\n";
    }
    else {
        cout << "Error: " << firstry.c_str() << " does not exist.  Please enter new file name\n";

        exit(1);
    }
}



//gabriel's version
void readICs_gen(string &firstry,  int &_Ntable3,Particle<2> *&_p,double const& efcheck, int & numpart)
{
    string filename;

    filename = ifolder+firstry;
    //nameenter:
    FILE * myfile = fopen (filename.c_str(),"r");

    if ((!firstry.empty())&&(myfile!= NULL)) {
        fscanf(myfile,"%i \n",&_Ntable3);
        _p= new Particle<2>[_Ntable3];
        int i=0;

        while (i<_Ntable3 ) {
            fscanf(myfile,"%lf    %lf    %lf\n",&_p[i].r.x[0],&_p[i].r.x[1],&_p[i].sigmaweight);
            _p[i].u.x[0]=0;
            _p[i].u.x[1]=0;
            _p[i].eta_sigma  = 1;
            _p[i].Bulk = 0;
            _p[i].Freeze=0;
            ++i;

        }


        fclose(myfile);
        cout << filename.c_str() << ": Input sucessful!\n";
    }
    else {
        cout << "Error: " << firstry.c_str() << " does not exist.  Please enter new file name\n";

        exit(1);
    }
}

//event by event
void readICs_ebe(string &firstry,  int &_Ntable3,Particle<2> *&_p,double factor,double const& efcheck, int & numpart)
{

    string filename;
    filename = ifolder+firstry;
    //nameenter:
    FILE *  myfile = fopen (filename.c_str(),"r");

    cout << "input from: " << filename.c_str() << endl;

    double *xsub,*ysub,*esub;
    if ((!firstry.empty())&&(myfile!= NULL)) {
        double stepx,stepy;
        fscanf(myfile,"%i   %lf    %lf %*f %*i %*i\n",&_Ntable3,&stepx,&stepy);
        // cout << "Initial Grid Size=" << _Ntable3 << endl;


        int i=0;



        xsub=new double[_Ntable3];
        ysub=new double[_Ntable3];
        esub=new double[_Ntable3];


        for(int j=0; j<_Ntable3; j++) {
            fscanf(myfile,"%lf    %lf    %lf \n",&xsub[i],&ysub[i],&esub[i]);


            esub[i]/=0.1973;
            // if (esub[i]>0.15) ++i;
            ++i;

        }


        _Ntable3=i;
        _p= new Particle<2>[_Ntable3];

        //  cout << "After e-cutoff=" << _Ntable3 << endl;


        int kk=_Ntable3;
        numpart=0;


        for(int j=0; j<_Ntable3; j++) {
            _p[j].r.x[0]=xsub[j];
            _p[j].r.x[1]=ysub[j];
            _p[j].e_sub=factor*esub[j];
            _p[j].u.x[0]=0;
            _p[j].u.x[1]=0;
            _p[j].eta_sigma  = 1;
            _p[j].sigmaweight=stepx*stepy;
            _p[j].Bulk = 0;



            if (_p[j].e_sub>efcheck)
            {
                _p[j].Freeze=0;
            }
            else
            {
                _p[j].Freeze=4;
                --kk;
                ++numpart;
            }
        }

        cout << "After freezeout=" << _Ntable3-numpart << endl;

        delete [] xsub;
        delete [] ysub;
        delete [] esub;




        fclose(myfile);
    }
    else {
        cout << "Error: " << firstry.c_str() << " does not exist.  Please enter new file name\n";

        exit(1);
    }

}


//trento
void readICs_tnt(string &firstry,  int &_Ntable3,Particle<2> *&_p,double factor,double const& sfcheck, int & numpart, eos EOS)
{

    string filename;
    filename = ifolder+firstry;
    //nameenter:
    ifstream input(filename.c_str());
    if (!input.is_open())
    {
        cout << "Can't open " << filename << endl;
        exit(1);
    }

    string line;
    vector<double> xsub,ysub,esub;

    getline(input,line);
    std::vector<std::string> gx = split(line, ' ');

    double stepx,stepy;
    stringstream s;
    s << gx[1];
    s >> stepx;


    stringstream s1;
    s1 << gx[2];
    s1 >> stepy;

    cout << "dx=dy=" << stepx << " " << stepy << endl;


    while (getline(input,line)) {
        std::vector<double> y (3,0) ;

        std::vector<std::string> x = split(line, ' ');


        for(int j=0; j<3; j++)
        {
            stringstream ss;
            ss << x[j];
            ss >> y[j];
//cout << "CHECK(" << __LINE__ << "): " << j << "   " << x[j] << "   " << y[j] << endl;
        }

        if ((factor*y[2])>0.001) {
            xsub.push_back(y[0]);
            ysub.push_back(y[1]);
            esub.push_back(y[2]);
        }

    }
    input.close();


    _Ntable3=xsub.size();
    _p= new Particle<2>[_Ntable3];

    cout << "After e-cutoff=" << _Ntable3 << endl;


    int kk=_Ntable3;
    numpart=0;



    for(int j=0; j<_Ntable3; j++) {
        _p[j].r.x[0]=xsub[j];
        _p[j].r.x[1]=ysub[j];
        // _p[j].e_sub=EOS.e_out(factor*esub[j]);
        _p[j].s_an=factor*esub[j];
        _p[j].u.x[0]=0;
        _p[j].u.x[1]=0;
        _p[j].eta_sigma  = 1;
        _p[j].sigmaweight=stepx*stepy;
        _p[j].Bulk = 0;



        if (_p[j].s_an>sfcheck)
        {
            _p[j].Freeze=0;
        }
        else
        {
            _p[j].Freeze=4;
            --kk;
            ++numpart;
        }
    }

    cout << "After freezeout=" << _Ntable3-numpart << endl;



}






//======================================================================
//iccing
void readICs_iccing( string &firstry, int &_Ntable3, Particle<2> *&_p,
					 double factor, double const & efcheck,
					 int & numpart, eos EOS)
{

	const double hbarC = 0.19733;
	cout << "Reading in ICCING initial conditions!" << endl;

    string filename;
    filename = ifolder+firstry;

	cout << "Initial conditions file: " << filename << endl;

    //nameenter:
    ifstream input(filename.c_str());
    if (!input.is_open())
    {
        cout << "Can't open " << filename << endl;
        exit(1);
    }

    string line;
    vector<double> xsub,ysub,esub,rBsub,rSsub,rQsub;

    getline(input,line);
    std::vector<std::string> gx = split(line, ' ');

    double stepx,stepy;
    stringstream s;
    s << gx[1];
    s >> stepx;


    stringstream s1;
    s1 << gx[2];
    s1 >> stepy;

    cout << "dx=dy=" << stepx << " " << stepy << endl;


    while (getline(input,line)) {
        std::vector<double> y (6,0) ;

        std::vector<std::string> x = split(line, ' ');


        for(int j=0; j<6; j++)
        {
            stringstream ss;
            ss << x[j];
            ss >> y[j];
//cout << "CHECK(" << __LINE__ << "): " << j << "   " << x[j] << "   " << y[j] << endl;
        }


		//==============================================================
		// Add new check here to enforce freeze-out criterion before
		// setting particle list size!!!

		// do not scale by factor!!!
        //if ((factor*y[2])>0.01)
        if (y[2]>0.01)
        //if (y[2]>max(0.01,0.5*efcheck*hbarC)) //N.B. leave wiggle room around FO temp
		{
			// check if this file has finite energy density and zero charge densities
			// (seems to cause problems in EOS)
			const double eps_local = 1e-6;
			if (   y[2] >= eps_local && y[3] < eps_local
				&& y[4] < eps_local  && y[5] < eps_local)
			{
				y[3] = eps_local;
				y[4] = eps_local;
				y[5] = eps_local;
			}
            xsub.push_back(y[0]);
            ysub.push_back(y[1]);
            esub.push_back(y[2]);
            rBsub.push_back(y[3]);
            rSsub.push_back(y[4]);
            rQsub.push_back(y[5]);
        }

    }

    input.close();

    _Ntable3=xsub.size();
    _p= new Particle<2>[_Ntable3];

    cout << "After e-cutoff and freeze-out: size = " << _Ntable3 << endl;


    int kk=_Ntable3;
    numpart=0;	//number of frozen out particles

    for(int j=0; j<_Ntable3; j++)
	{
        _p[j].r.x[0]=xsub[j];
        _p[j].r.x[1]=ysub[j];
        // _p[j].e_sub=EOS.e_out(factor*esub[j]);
        _p[j].e_sub=esub[j]/hbarC;        // not s_an!!  convert to 1/fm^4 and do not rescale by factor!
        _p[j].u.x[0]=0;
        _p[j].u.x[1]=0;
        _p[j].eta_sigma = 1;
        _p[j].sigmaweight=stepx*stepy;
		_p[j].rho_weight = stepx*stepy;
        _p[j].Bulk = 0;
        _p[j].B=rBsub[j]*stepx*stepy;			// confirm with Jaki
        _p[j].S=rSsub[j]*stepx*stepy;			// confirm with Jaki
        _p[j].Q=rQsub[j]*stepx*stepy;			// confirm with Jaki
        _p[j].rhoB_an=rBsub[j];					// confirm with Jaki
        _p[j].rhoS_an=rSsub[j];					// confirm with Jaki
        _p[j].rhoQ_an=rQsub[j];					// confirm with Jaki
		_p[j].transverse_area = stepx*stepy;

		//if (j==0)
		//cout << "readICs_iccing(" << __LINE__ << "): "
		cout << "SPH particles: "
			<< j << "   "
			<< _p[j].r.x[0] << "   " << _p[j].r.x[1] << "   "
			<< _p[j].e_sub << "   " << _p[j].rhoB_an << "   "
			<< _p[j].rhoS_an << "   " << _p[j].rhoQ_an << "   "
			<< _p[j].sigmaweight << endl;

		// make educated initial guess here for this particle's (T, mu_i) coordinates
		// (improve this in the future)
		_p[j].particle_T   = 1150.0/197.33;	// rootfinder seems to work better going downhill than "uphill"
		_p[j].particle_muB = 0.0;
		_p[j].particle_muS = 0.0;
		_p[j].particle_muQ = 0.0;

        if (_p[j].e_sub>efcheck)	// impose freeze-out check for e, not s
        {
            _p[j].Freeze=0;
        }
        else
        {
            _p[j].Freeze=4;
            --kk;
            ++numpart;
        }
    }

    cout << "After freezeout (redundant): size = " << _Ntable3-numpart << endl;



}


void readICs_iccing( string &firstry, int &_Ntable3, Particle<3> *&_p,
                     double factor, double const & efcheck, int & numpart, eos EOS)
{
}


//event by event for giorgio
void readICs_gebe(string &firstry,  int &_Ntable3,Particle<2> *&_p,double factor,double const& efcheck, int & numpart)
{

    string filename;
    filename = ifolder+firstry;
    //nameenter:
    ifstream input(filename.c_str());

    cout << "input from: " << filename.c_str() << endl;


    if (!input.is_open())
    {
        cout << "Can't open " << filename << endl;
        exit(1);
    }



    string line;
    vector<double> xsub,ysub,esub;

    double gd2;
    cout << "hcor=" <<  hcor << endl;
    if (hcor==1) {
        getline(input,line);
        cout << "here" << endl;
        gd2=0.06*0.06;
    }
    else {
        gd2=0.08*0.08;
    }


    while (getline(input,line)) {
        std::vector<double> y (3,0) ;

        std::vector<std::string> x = split(line, ' ');


        for(int j=0; j<3; j++)
        {
            stringstream ss;
            ss << x[j];
            ss >> y[j];
        }

        if (y[2]>0.01) {
            xsub.push_back(y[0]);
            ysub.push_back(y[1]);
            esub.push_back(y[2]);
        }

    }
    input.close();


    _Ntable3=xsub.size();
    _p= new Particle<2>[_Ntable3];

    cout << "After e-cutoff=" << _Ntable3 << endl;


    int kk=_Ntable3;
    numpart=0;



    for(int j=0; j<_Ntable3; j++) {
        _p[j].r.x[0]=xsub[j];
        _p[j].r.x[1]=ysub[j];
        _p[j].e_sub=esub[j];
        _p[j].u.x[0]=0;
        _p[j].u.x[1]=0;
        _p[j].eta_sigma  = 1;
        _p[j].sigmaweight=gd2;
        _p[j].Bulk = 0;



        if (_p[j].e_sub>efcheck)
        {
            _p[j].Freeze=0;
        }
        else
        {
            _p[j].Freeze=4;
            --kk;
            ++numpart;
        }
    }

    cout << "After freezeout=" << _Ntable3-numpart << endl;

}




//event by event for nexus with flow
void readICs_nebe(string &firstry,  int &_Ntable3,Particle<2> *&_p,double factor,double const& efcheck, int & numpart)
{

    string filename;
    filename = ifolder+firstry;
    //nameenter:
    FILE *  myfile = fopen (filename.c_str(),"r");

    cout << "input from: " << filename.c_str() << endl;
    _Ntable3=50000;
    double *xsub,*ysub,*esub,*vxsub,*vysub,*vzsub;
    if ((!firstry.empty())&&(myfile!= NULL)) {


        int i=0;



        xsub=new double[_Ntable3];
        ysub=new double[_Ntable3];
        esub=new double[_Ntable3];
        vxsub=new double[_Ntable3];
        vysub=new double[_Ntable3];
        vzsub=new double[_Ntable3];


        while (fscanf(myfile,"%lf    %lf  %lf %lf    %lf  %lf  \n",&xsub[i],&ysub[i],&esub[i],&vxsub[i],&vysub[i],&vzsub[i])!=EOF) {


            esub[i]/=0.1973;
            if (factor*esub[i]>0.01) ++i;
        }



        _Ntable3=i;
        _p= new Particle<2>[_Ntable3];

        // cout << "After e-cutoff=" << _Ntable3 << endl;

        int kk=_Ntable3;
        numpart=0;

        for(int j=0; j<_Ntable3; j++) {
            _p[j].r.x[0]=xsub[j];
            _p[j].r.x[1]=ysub[j];
            _p[j].e_sub=factor*esub[j];
            _p[j].v.x[0]=vxsub[j];
            _p[j].v.x[1]=vysub[j];
            _p[j].gamma=1/sqrt(1-Norm2(_p[j].v));
            _p[j].u.x[0]=_p[j].gamma*vxsub[j];
            _p[j].u.x[1]=_p[j].gamma*vysub[j];

            _p[j].eta_sigma  = 1;
            _p[j].sigmaweight=0.08*0.08;
            _p[j].Bulk = 0;



            if (_p[j].e_sub>efcheck)
            {
                _p[j].Freeze=0;
            }
            else
            {
                _p[j].Freeze=4;
                --kk;
                ++numpart;
            }
        }
        cout << "After freezeout=" << _Ntable3-numpart << endl;

        delete [] xsub;
        delete [] ysub;
        delete [] esub;

        fclose(myfile);
    }
    else {
        cout << "Error: " << firstry.c_str() << " does not exist.  Please enter new file name\n";

        exit(1);
    }
}

//event by event for ipglasma with flow
double readICs_gl(string &firstry,  int &_Ntable3,Particle<2> *&_p,double factor,double const& efcheck, int & numpart)
{

    string filename;
    filename = ifolder+firstry;


    ifstream input(filename.c_str());

    cout << "input from: " << filename.c_str() << endl;

    vector< double > rx,ry,e,ux,uy;

    if (!input.is_open())
    {
        cout << "Can't open " << filename << endl;
        exit(1);
    }




    double gd2;
    if (hcor==2) gd2=0.17*0.17/4.;
    else if (hcor==4) gd2=0.17*0.17/16.;
    else if (hcor==6) {
        gd2=0.0283333*0.0283333;

        string line;
        getline(input,line);
        while (getline(input,line)) {

            std::vector<double> y (7,0) ;
            std::vector<std::string> x = split(line, ' ');

            if (x.size()<7) break;

            for(int j=0; j<7; j++)
            {
                stringstream s;
                s << x[j];
                s >> y[j];
            }

            double esub= factor*y[3];
            if (esub>0.05)
            {
                e.push_back(esub);
                rx.push_back(y[1]);
                ry.push_back(y[2]);
                ux.push_back(y[4]*y[5]);
                uy.push_back(y[4]*y[6]);

            }


        }


    }
    else if (hcor==7) {
        gd2=0.0283333*0.0283333*4;

        string line;
        getline(input,line);
        int skip=0,run=0;
        while (getline(input,line)) {

            std::vector<double> y (7,0) ;
            std::vector<std::string> x = split(line, ' ');

            run++;
            if (run==1200) skip=0;
            if (run==2401) {
                skip=0;
                run=1;
            }
            if (skip==1) {
                skip=0;
                continue;
            }
            else skip=1;

            if (run>1200) continue;


            if (x.size()<7) break;

            for(int j=0; j<7; j++)
            {
                stringstream s;
                s << x[j];
                s >> y[j];
            }
//       cout << run << " " << skip << " " <<  x[1] << " " << x[2] << endl;
//       if (run==4) getchar();
//       if (run==5) getchar();

            double esub= factor*y[3];
            if (esub>0.05)
            {
                e.push_back(esub);
                rx.push_back(y[1]);
                ry.push_back(y[2]);
                ux.push_back(y[4]*y[5]);
                uy.push_back(y[4]*y[6]);

            }


        }
    }
    else  gd2=0.17*0.17;



    if (hcor<6) {

        string line;
        getline(input,line);


        while (getline(input,line)) {


            std::vector<double> y (6,0) ;
            std::vector<std::string> x = split(line, ' ');

            if (x.size()<5) break;

            for(int j=0; j<5; j++)
            {
                stringstream s;
                s << x[j];
                s >> y[j];
            }


            double esub= factor*y[2]/0.1973;
//       if (esub>0.1)
//       {



            e.push_back(esub);
            rx.push_back(y[0]);
            ry.push_back(y[1]);
            ux.push_back(y[3]);
            uy.push_back(y[4]);
//       }

        }
    }


    input.close();



    _Ntable3=e.size();
    _p= new Particle<2>[_Ntable3];


    int kk=_Ntable3;
    numpart=0;

    for(int j=0; j<_Ntable3; j++) {
        _p[j].r.x[0]=rx[j];
        _p[j].r.x[1]=ry[j];
        _p[j].e_sub=e[j];
        _p[j].u.x[0]=ux[j];  //these are actually already u not v
        _p[j].u.x[1]=uy[j];

        _p[j].eta_sigma  = 1;
        _p[j].sigmaweight=gd2;
        _p[j].Bulk = 0;



        if (_p[j].e_sub>efcheck)
        {
            _p[j].Freeze=0;
        }
        else
        {
            _p[j].Freeze=4;
            --kk;
            ++numpart;
        }
    }
    cout << "After freezeout=" << _Ntable3-numpart << endl;



    return gd2;

}

//event by event for ipglasma with flow
void readICs_glno(string &firstry,  int &_Ntable3,Particle<2> *&_p,double factor,double const& efcheck, int & numpart)
{

    string filename;
    filename = ifolder+firstry;


    ifstream input(filename.c_str());

    cout << "input from: " << filename.c_str() << endl;

    vector< double > rx,ry,e;

    if (!input.is_open())
    {
        cout << "Can't open " << filename << endl;
        exit(1);
    }




    double gd2;
    if (hcor==2) gd2=0.17*0.17/4.;
    else if (hcor==4) gd2=0.17*0.17/16.;
    else if (hcor==6) {
        gd2=0.0283333*0.0283333;

        string line;
        getline(input,line);
        while (getline(input,line)) {

            std::vector<double> y (7,0) ;
            std::vector<std::string> x = split(line, ' ');

            if (x.size()<7) break;

            for(int j=0; j<7; j++)
            {
                stringstream s;
                s << x[j];
                s >> y[j];
            }

            double esub= factor*y[3];
            if (esub>0.05)
            {
                e.push_back(esub);
                rx.push_back(y[1]);
                ry.push_back(y[2]);

            }


        }


    }
    else if (hcor==7) {
        gd2=0.0283333*0.0283333*4;

        string line;
        getline(input,line);
        int skip=0,run=0;
        while (getline(input,line)) {

            std::vector<double> y (7,0) ;
            std::vector<std::string> x = split(line, ' ');

            run++;
            if (run==1200) skip=0;
            if (run==2401) {
                skip=0;
                run=1;
            }
            if (skip==1) {
                skip=0;
                continue;
            }
            else skip=1;

            if (run>1200) continue;


            if (x.size()<7) break;

            for(int j=0; j<7; j++)
            {
                stringstream s;
                s << x[j];
                s >> y[j];
            }
//       cout << run << " " << skip << " " <<  x[1] << " " << x[2] << endl;
//       if (run==4) getchar();
//       if (run==5) getchar();

            double esub= factor*y[3];
            if (esub>0.05)
            {
                e.push_back(esub);
                rx.push_back(y[1]);
                ry.push_back(y[2]);

            }


        }
    }
    else  gd2=0.17*0.17;



    if (hcor!=6) {
        string line;

        getline(input,line);
        while (getline(input,line)) {

            std::vector<double> y (5,0) ;
            std::vector<std::string> x = split(line, ' ');

            if (x.size()<5) break;

            for(int j=0; j<5; j++)
            {
                stringstream s;
                s << x[j];
                s >> y[j];
            }

            double esub= factor*y[2]/0.1973;
            if (esub>0.1)
            {

                e.push_back(esub);
                rx.push_back(y[0]);
                ry.push_back(y[1]);
            }

        }
    }


    input.close();



    _Ntable3=e.size();
    _p= new Particle<2>[_Ntable3];


    int kk=_Ntable3;
    numpart=0;

    for(int j=0; j<_Ntable3; j++) {
        _p[j].r.x[0]=rx[j];
        _p[j].r.x[1]=ry[j];
        _p[j].e_sub=e[j];
        _p[j].u.x[0]=0;  //these are actually already u not v
        _p[j].u.x[1]=0;

        _p[j].eta_sigma  = 1;
        _p[j].sigmaweight=gd2;
        _p[j].Bulk = 0;



        if (_p[j].e_sub>efcheck)
        {
            _p[j].Freeze=0;
        }
        else
        {
            _p[j].Freeze=4;
            --kk;
            ++numpart;
        }
    }
    cout << "After freezeout=" << _Ntable3-numpart << endl;
}

double readICs_gl(string &firstry,  int &_Ntable3,Particle<3> *&_p,double factor,double const& efcheck, int & numpart) {

    return factor;
}

void readICs_glno(string &firstry,  int &_Ntable3,Particle<3> *&_p,double factor,double const& efcheck, int & numpart) {
}

void readICs_gebe(string &firstry,  int &_Ntable3,Particle<3> *&_p,double factor,double const& efcheck, int & numpart)
{
}

void readICs_nebe(string &firstry,  int &_Ntable3,Particle<3> *&_p,double factor,double const& efcheck, int & numpart)
{
}

void readICs_tnt(string &firstry,  int &_Ntable3,Particle<3> *&_p,double factor,double const& efcheck, int & numpart, eos EOS) {
}

//with initial flow
void readICs_genUU(string &firstry,  int &_Ntable3,Particle<2> *&_p,double const& efcheck, int & numpart)
{
    FILE * myfile;
    string filename;
    int i;

    filename = ifolder+firstry;
    //nameenter:
    myfile = fopen (filename.c_str(),"r");
    if ((!firstry.empty())&&(myfile!= NULL)) {
        fscanf(myfile,"%i \n",&_Ntable3);
        _p= new Particle<2>[_Ntable3];
        i=0;

        while (i<_Ntable3 ) {
            fscanf(myfile,"%lf    %lf    %lf  %lf    %lf\n",&_p[i].r.x[0],&_p[i].r.x[1],&_p[i].sigmaweight,&_p[i].u.x[0],&_p[i].u.x[1]);
            _p[i].eta_sigma  = 1.;
            _p[i].Bulk = 0.;

            if (_p[i].e_sub>efcheck)
            {
                _p[i].Freeze=0;
            }
            else
            {
                _p[i].Freeze=4;
                numpart++;
            }

            i++;
        }
        fclose(myfile);
        cout << filename.c_str() << ": Input sucessful!\n";
    }
    else {
        cout << "Error: " << firstry.c_str() << " does not exist.  Please enter new file name\n";

        exit(1);
    }
}

//3 dimensions

void readICs_gen(string &firstry,  string &firstry2, int &_Ntable3, Particle< 3 > *&_p,double const& efcheck, int & numpart)
{
    FILE * myfile;
    FILE * myfile2;
    string filename, filename2;
    int i;

    filename = ifolder+firstry;
    //nameenter:
    myfile = fopen (filename.c_str(),"r");
    if ((!firstry.empty())&&(myfile!= NULL)) {
        fscanf(myfile,"%i \n",&_Ntable3);
        _p= new Particle<3>[_Ntable3];
        i=0;
        while (i<_Ntable3 ) {
            fscanf(myfile,"%lf    %lf    %lf   %lf\n",&_p[i].r.x[0],&_p[i].r.x[1],&_p[i].r.x[2],&_p[i].sigmaweight);
            _p[i].eta_sigma  = 1.;
            _p[i].Bulk = 0.;
            i++;
        }

        fclose(myfile);
        cout << filename.c_str() << ": Input sucessful!\n";
    }
    else {
        cout << "Error: " << firstry.c_str() << " does not exist.  Please enter new file name\n";

        exit(1);
    }



    filename2 = ifolder+firstry2;
    //nameenter2:
    myfile2 = fopen (filename2.c_str(),"r");
    if ((!firstry2.empty())&&(myfile2!= NULL)) {

        i=0;
        while (i<_Ntable3 ) {
            fscanf(myfile2,"%*f  %*f  %*f    %lf    %lf   %lf\n",&_p[i].v.x[0],&_p[i].v.x[1],&_p[i].v.x[2]);
//                   _p[i].gammacalc();
//                   _p[i].ucalc();
            i++;
        }

        fclose(myfile2);
        cout << filename2.c_str() << ": Input sucessful!\n";
    }
    else {
        cout << "Error: " << firstry2.c_str() << " does not exist.  Please enter new file name\n";

        exit(1);
    }
}

void readICs_gen(string &firstry,  int &_Ntable3, Particle<3> *&_p,double const& efcheck, int & numpart)
{
    FILE * myfile;
    string filename;
    int i;


    filename = ifolder+firstry;
    //nameenter:
    myfile = fopen (filename.c_str(),"r");
    if ((!firstry.empty())&&(myfile!= NULL)) {
        fscanf(myfile,"%i \n",&_Ntable3);
        _p= new Particle<3>[_Ntable3];
        i=0;
        while (i<_Ntable3 ) {
            fscanf(myfile,"%lf    %lf    %lf   %lf\n",&_p[i].r.x[0],&_p[i].r.x[1],&_p[i].r.x[2],&_p[i].sigmaweight);
            _p[i].eta_sigma  = 1.;
            _p[i].u.x[0]=0;
            _p[i].u.x[1]=0;
            _p[i].u.x[2]=0;
            _p[i].Bulk = 0.;
            i++;
        }

        fclose(myfile);
        cout << filename.c_str() << ": Input sucessful!\n";
    }
    else {
        cout << "Error: " << firstry.c_str() << " does not exist.  Please enter new file name\n";

        exit(1);
    }



}

//with initial flow
void readICs_genUU(string &firstry,  int &_Ntable3,Particle<3> *&_p,double const& efcheck, int & numpart)
{
    FILE * myfile;
    string filename;
    int i;


    filename = ifolder+firstry;
    //nameenter:
    myfile = fopen (filename.c_str(),"r");
    if ((!firstry.empty())&&(myfile!= NULL)) {
        fscanf(myfile,"%i \n",&_Ntable3);
        _p= new Particle<3>[_Ntable3];
        i=0;

        while (i<_Ntable3 ) {
            fscanf(myfile,"%lf    %lf    %lf  %lf    %lf %lf    %lf\n",&_p[i].r.x[0],&_p[i].r.x[1],&_p[i].r.x[2],&_p[i].sigmaweight,&_p[i].u.x[0],&_p[i].u.x[1],&_p[i].u.x[2]);
            _p[i].eta_sigma  = 1.;
            _p[i].Bulk = 0.;


            i++;
        }
        fclose(myfile);
        cout << filename.c_str() << ": Input sucessful!\n";
    }
    else {
        cout << "Error: " << firstry.c_str() << " does not exist.  Please enter new file name\n";

        exit(1);
    }
}


void open_dir(string eventfolder, string *&filelist, int &count)
{

    string dotslash ("./");
    string dash ("/");
    string dir_in=dotslash+ifolder+("/")+eventfolder+("/");
    string dat("dat");
    string d1 ("dens1.dat");
    string d2 ("dens2.dat");
    DIR *pdir;
    count=0;
    struct dirent *pent;
    pdir=opendir(dir_in.c_str());
    if (!pdir)
    {
        cout<<"Directory doesnot exist";
    }
    int errno4=0; //errno.h
    while ((pent=readdir(pdir)))
    {
        string file_open=pent->d_name;
        int j=file_open.length();
        if (j>4)
        {
            string datcheck=file_open.substr((j-3),3);
            if (dat==datcheck&&file_open!=d1&&file_open!=d2) count++;
        }


    }
    closedir(pdir);

    filelist=new string[count];
    pdir=opendir(dir_in.c_str());
    int i=0;
    while ((pent=readdir(pdir)))
    {
        string file_open=pent->d_name;
        int j=file_open.length();

        if (j>4)
        {
            string datcheck=file_open.substr((j-3),3);
            if (dat==datcheck&&file_open!=d1&&file_open!=d2)
            {
                filelist[i]= eventfolder+dash+file_open;
                i++;
            }
        }
    }
    if (errno4)
    {
        cout<<"Error while accessing directory";
    }
    closedir(pdir);

}

////event by event
void readICs_ebe(string &firstry,  int &_Ntable3,Particle<3> *&_p,double factor, double const& efcheck, int & numpart)
{
    FILE * myfile;
    string filename;
    int i;
    double stepx,stepy,stepeta;


    filename = ifolder+firstry;
    //nameenter:
    myfile = fopen (filename.c_str(),"r");
    if ((!firstry.empty())&&(myfile!= NULL)) {
        fscanf(myfile,"%i %lf    %lf   %lf\n",&_Ntable3,&stepx,&stepy,&stepeta);
        _p= new Particle<3>[_Ntable3];
        i=0;


        while (i<_Ntable3 ) {
            fscanf(myfile,"%lf    %lf    %lf   %lf\n",&_p[i].r.x[0],&_p[i].r.x[1],&_p[i].r.x[2],&_p[i].e_sub);
            _p[i].u.x[0]=0.;
            _p[i].u.x[1]=0.;
            _p[i].u.x[2]=0.;
            _p[i].e_sub/=0.1973;
            _p[i].eta_sigma  = 1.;
            _p[i].sigmaweight  =stepx*stepy*stepeta;
            _p[i].Bulk = 0.;

            i++;
        }





        fclose(myfile);
        cout << filename.c_str() << ": Input sucessful!\n";
    }
    else {
        cout << "Error: " << firstry.c_str() << " does not exist.  Please enter new file name\n";

        exit(1);
    }
}


string convertInt(int number)
{
    stringstream ss;//create a stringstream
    ss << number;//add number to the stream
    return ss.str();//return a string with the contents of the stream
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if (!item.empty()) elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

#endif
