#ifndef _OUTPUT_H_
#define _OUTPUT_H_

#include<fstream>
#include<stdio.h>
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

using namespace std;

template <int D>
class Output {
private:
    ofstream ESC,FO,SC;
    string ESCname,FOname;
    string Scons;
    int countEP,countGC;
    string convertInt(int number);
    string cfolder;
    typedef struct stat Stat;

public:
    Output<D>(LinkList<D> &linklist);
    void cleardir(LinkList<D> &linklist);
    void conservationstart(LinkList<D> &linklist);
    void conservation(LinkList<D> &linklist);
//    void sconservation(LinkList<D> &linklist);
    void eprofile(LinkList<D> &linklist);
    void eccout(LinkList<D> &linklist);
    double ecc(LinkList<D> &linklist, double  & psi,double & rout, int m,int n);
    void sveprofile(LinkList<D> &linklist);
	void print_physical_quantities( LinkList<D> &linklist );
	void compute_physical_quantities( LinkList<D> &linklist,
		Vector<int,D> r0, double & temperature, double & baryon_chemical_potential,
		double & strange_chemical_potential, double & electric_chemical_potential,
		double & energy_density, double & baryon_density, double & strange_density,
		double & electric_density, double & entropy_density );
    void bsqsveprofile(LinkList<D> &linklist);
    void gubcheckux(LinkList<D> &linklist);
    void gubcheckuy(LinkList<D> &linklist);
    void do_mkdir(const char *path, int mode);
    void averages(LinkList<D> &linklist);
    void FOstart(LinkList<D> &linklist);
    void FOprint(LinkList<D> &linklist);
    void vFOprint(LinkList<D> &linklist);
    void svFOprint(LinkList<D> &linklist);
    void bsqsvFOprint(LinkList<D> &linklist);
    void sconstart(LinkList<D> &linklist);
    void SCprint(LinkList<D> &linklist);
};

template <int D>
Output<D>::Output(LinkList<D> &linklist)
{
    if (linklist.fcount!=0)
    {
        cfolder=ofolder+linklist.ebe_folder+("/");
        if (linklist.cevent==0)
        {
            if (linklist.fcount!=0) do_mkdir(cfolder.c_str(), 0777);
            cleardir(linklist);
        }
        cout << "output to:" << cfolder << endl;

    }

    if (linklist.average==1)
    {
        cfolder=ofolder+linklist.ebe_folder+("/");
        do_mkdir(cfolder.c_str(), 0777);
    }


    countEP=0;
    // conservationstart( linklist);
    FOstart( linklist);
    if  (linklist.visc==1) sconstart(linklist);


}

template <int D>
void Output<D>::conservationstart(LinkList<D> &linklist)
{


    string esin;

    if (linklist.average==1)
    {
        ESCname=cfolder+"average_Econsv.dat";
    }
    else if (linklist.fcount==0&&linklist.average!=1)
    {
        if (linklist.visc==0)
        {
            esin="conserv";
        }
        else if  (linklist.visc==1)
        {
            esin="bvconserv";
        }
        else if  (linklist.visc==3)
        {
            esin="sbvconserv";
        }


        ESCname=ofolder+esin+"_weight.dat";
    }
    else
    {

        if (linklist.visc==0)
        {
            esin="conserv";
        }
        else if  (linklist.visc==1)
        {
            esin="bvconserv";
        }
        else if  (linklist.visc==3)
        {
            esin="sbvconserv";
        }

        string under ("_ev");
        string even;
        even=convertInt(linklist.fnum);
        ESCname=cfolder+esin+under+even+".dat";
    }




    ESC.open(ESCname.c_str());
    if (!ESC.is_open())
    {

        cout << "Error: cannot open E_S_conserv.dat file!" << endl;
        exit(1);
    }

    ESC.close();

    countEP=0;
    countGC=0;

}

template <int D>
void Output<D>::FOstart(LinkList<D> &linklist)
{


    string FOin;

    if (linklist.average==1)
    {
        FOname=cfolder+"average_freezeout.dat";
    }
    else if (linklist.fcount==0&&linklist.average!=1)
    {
        if (linklist.visc==0)
        {
            FOin="freezeout";
        }
        else if  (linklist.visc==1)
        {
            FOin="bvfreezeout";
        }
        else if  (linklist.visc==3)
        {
            FOin="sbvfreezeout";
        }




        FOname=ofolder+FOin+"_weight.dat";
    }
    else
    {

        string under ("_ev");
        string even;
        even=convertInt(linklist.fnum);
        if (linklist.visc==0)
        {
            FOin="freezeout";
        }
        else if  (linklist.visc==1)
        {
            FOin="bvfreezeout";
        }
        else if  (linklist.visc==3)
        {
            FOin="sbvfreezeout";
        }


        FOname=cfolder+FOin+under+even+".dat";
    }




    FO.open(FOname.c_str());
    if (!FO.is_open())
    {

        cout << "Error: cannot open E_S_conserv.dat file!" << endl;
        exit(1);
    }
//    if  (linklist.visc==3) {//
//        //double therm=1./(linklist.sFO*pow(freezeoutT,3));
//        FO << linklist.sFO << endl;
//
//    }
//    else  FO << linklist.sFO << endl;

    FO.close();


}

template <int D>
void Output<D>::sconstart(LinkList<D> &linklist)
{
    Scons=cfolder+"scons.dat";

    if (linklist.fcount!=0&&linklist.average!=1)
    {
        string even;
        even=convertInt(linklist.fnum);
        Scons=cfolder+"scons"+even+".dat";
    }

    SC.open(Scons.c_str());
    if (!SC.is_open())
    {

        cout << "Error: cannot open scons.dat file!" << endl;
        exit(1);
    }

    SC.close();
}

template <int D>
void Output<D>::gubcheckux(LinkList<D> &linklist)
{
    ofstream GSC;

    countGC+=1;
    string gcin ("gubcheck_ux");
    string dat (".dat");
    string sCEP;
    sCEP=convertInt(countEP);

    string gcname=ofolder+gcin+sCEP+dat;


    GSC.open(gcname.c_str());
    if (!GSC.is_open())
    {

        cout << "Error: cannot open gubcheck" << sCEP <<".dat file!" << endl;
        exit(1);
    }
    else
    {
        //GSC << linklist.t << endl;
        for (int i=0; i<linklist.n(); i++)
        {
            if (abs(Norm(linklist._p[i].r ))<4.)
            {
                GSC << linklist._p[i].r  <<  " " << linklist._p[i].u.x[0] << endl;
            }
            //GSC << linklist._p[i].r  <<  " " << linklist._p[i].u.x[0]/linklist.ux(linklist.t,linklist._p[i].r.x[0],linklist._p[i].r.x[1]) << endl;

        }
    }



    GSC.close();

}

template <int D>
void Output<D>::gubcheckuy(LinkList<D> &linklist)
{
    ofstream GSC;

    countGC+=1;
    string gcin ("gubcheck_uy");
    string dat (".dat");
    string sCEP;
    sCEP=convertInt(countEP);
    string gcname=ofolder+gcin+sCEP+dat;





    GSC.open(gcname.c_str());
    if (!GSC.is_open())
    {

        cout << "Error: cannot open gubcheck" << sCEP <<".dat file!" << endl;
        exit(1);
    }
    else
    {
        //GSC << linklist.t << endl;
        for (int i=0; i<linklist.n(); i++)
        {
            //    GSC << linklist._p[i].r  <<  " " << linklist._p[i].u.x[1]/linklist.uy(linklist.t,linklist._p[i].r.x[0],linklist._p[i].r.x[1]) << endl;
            if (abs(Norm(linklist._p[i].r ))<4.)
            {
                GSC << linklist._p[i].r  <<  " " << linklist._p[i].u.x[1] << endl;
            }
        }
    }



    GSC.close();

}

template <int D>
void Output<D>::eprofile(LinkList<D> &linklist)
{
    ofstream EPN;

    countEP+=1;
    string epin ("eprofile");
    string dat (".dat");
    string sCEP;
    sCEP=convertInt(countEP);
    string epname;

    if (linklist.average==1)
    {
        epname=cfolder+"average_eprof"+ sCEP + ".dat";
    }
    else if (linklist.fcount==0&&linklist.average!=1)
    {
        epname=ofolder+epin+sCEP+dat;
    }
    else
    {
        string under ("_ev");
        string even;
        even=convertInt(linklist.fnum);
        epname=cfolder+epin+sCEP+under+even+dat;
    }

    //cout << epname << endl;
    EPN.open(epname.c_str());
    if (!EPN.is_open())
    {

        cout << "Error: cannot open eprofile" << sCEP <<".dat file!" << endl;
        exit(1);
    }
    else
    {
        EPN << linklist.t << endl;


        for (int i=0; i<linklist.n(); i++)
        {
//            for (int j=0;j<linklist.D;j++)
//            {
//                 EPN << linklist._p[i].r.x[j] << " ";
//             }

//    EPN << linklist._p[i].r   << " " << linklist._p[i].EOSe()/linklist.eanal2(linklist.t,linklist._p[i].r.x[0],linklist._p[i].r.x[1]) <<  endl;
//    if (abs(Norm(linklist._p[i].r ))<4.)
//        {

            EPN << linklist._p[i].r   << " " << linklist._p[i].EOSe() << " " << linklist._p[i].EOST()*197.3 <<  endl;
//        EPN << linklist.t   << " "  << linklist._p[i].r   << " " << linklist._p[i].EOST()*0.1973 << " " <<  linklist._p[i].EOSp()   << " " <<  linklist._p[i].EOSe()   << " " <<  linklist._p[i].EOSs()   << " " <<  linklist._p[i].u   <<endl;
            //}
            //cout << linklist._p[i].EOSs() << " " << linklist._p[i].sigma << endl;
            //getchar();

        }
    }



    EPN.close();

}








template <int D>
void Output<D>::sveprofile(LinkList<D> &linklist)
{
    ofstream EPN;

    countEP+=1;
    string epin ("sveprofile");
    string dat (".dat");
    string sCEP;
    sCEP=convertInt(countEP);
    string epname;

    if (linklist.average==1)
    {
        epname=cfolder+"average_eprof"+ sCEP + ".dat";
    }
    else if (linklist.fcount==0&&linklist.average!=1)
    {
        epname=ofolder+epin+sCEP+dat;
    }
    else
    {
        string under ("_ev");
        string even;
        even=convertInt(linklist.fnum);
        epname=cfolder+epin+sCEP+under+even+dat;
    }

    //cout << epname << endl;
    EPN.open(epname.c_str());
    if (!EPN.is_open())
    {

        cout << "Error: cannot open eprofile" << sCEP <<".dat file!" << endl;
        exit(1);
    }
    else
    {
        EPN << linklist.t << endl;



        for (int i=0; i<linklist.n(); i++)
        {

//            for (int j=0;j<linklist.D;j++)
//            {
//                 EPN << linklist._p[i].r.x[j] << " ";
//             }

//    EPN << linklist._p[i].r   << " " << linklist._p[i].EOSe()/linklist.eanal2(linklist.t,linklist._p[i].r.x[0],linklist._p[i].r.x[1]) <<  endl;
            //if (abs(Norm(linklist._p[i].r ))<4.)
            //{






//        EPN << linklist._p[i].r   << " " << linklist._p[i].v  << " " <<  linklist._p[i].shv.x[0][0]*0.1973  << " " << linklist._p[i].shv.x[1][1]*0.1973  << " " << linklist._p[i].shv.x[1][2]*0.1973  << " " << linklist._p[i].shv.x[2][2]*0.1973  << " " << linklist.t*linklist.t*linklist._p[i].shv33*0.1973  << " " << linklist._p[i].EOSp()*0.1973  <<  " " <<linklist._p[i].EOSe()*0.1973 <<   " " << linklist._p[i].EOST()*0.1973 << endl;



            //EPN << linklist._p[i].EOSe()<< " " << linklist._p[i].EOSe() << " " << linklist._p[i].EOST()*197.3 << " " << linklist._p[i].v  << " " << linklist._p[i].sigmaweight   <<endl;

            EPN << linklist._p[i].r   << " " << linklist._p[i].EOSe() << " " << linklist._p[i].EOSp() << " " << linklist._p[i].EOST()*197.3 << " " << linklist._p[i].stauRelax << " " <<  linklist._p[i].sigmaweight << " " << linklist._p[i].bigtheta   << " " <<  sqrt(linklist._p[i].shv.x[0][0]*linklist._p[i].shv.x[0][0]-2*linklist._p[i].shv.x[0][1]*linklist._p[i].shv.x[0][1] -2*linklist._p[i].shv.x[0][2]*linklist._p[i].shv.x[0][2]    + linklist._p[i].shv.x[1][1]*linklist._p[i].shv.x[1][1]+ linklist._p[i].shv.x[2][2]*linklist._p[i].shv.x[2][2] +2* linklist._p[i].shv.x[1][2]*linklist._p[i].shv.x[1][2]+pow(linklist.t,4)*linklist._p[i].shv33*linklist._p[i].shv33) << " " << linklist._p[i].stauRelax/linklist.t * linklist._p[i].bigtheta << " " << linklist._p[i].u.x[0]/linklist._p[i].gamma << " " << linklist._p[i].u.x[1]/linklist._p[i].gamma << " " << linklist._p[i].gamma <<endl;

//        EPN << linklist._p[i].r   << " " << linklist._p[i].stauRelax/linklist.t * linklist._p[i].bigtheta  << " " <<  sqrt(linklist._p[i].shv.x[0][0]*linklist._p[i].shv.x[0][0]-2*linklist._p[i].shv.x[0][1]*linklist._p[i].shv.x[0][1] -2*linklist._p[i].shv.x[0][2]*linklist._p[i].shv.x[0][2]    + linklist._p[i].shv.x[1][1]*linklist._p[i].shv.x[1][1]+ linklist._p[i].shv.x[2][2]*linklist._p[i].shv.x[2][2] +2* linklist._p[i].shv.x[1][2]*linklist._p[i].shv.x[1][2]+pow(linklist.t,4)*linklist._p[i].shv33*linklist._p[i].shv33)/linklist._p[i].EOSp() << " " << linklist._p[i].EOST()*0.1973 << "  " <<  linklist._p[i].EOSe() << " " << linklist._p[i].u.x[0]/linklist._p[i].gamma << " " << linklist._p[i].u.x[1]/linklist._p[i].gamma <<  endl;


            //cout << linklist._p[i].EOSs() << " " << linklist._p[i].sigma << endl;
            //getchar();

        }
    }



    EPN.close();

}


///*
template <int D>
void Output<D>::print_physical_quantities( LinkList<D> &linklist )
{
	const double dx = 0.1, dy = 0.1;
	const double xmin = -15.0, ymin = -15.0;
	const double xmax = -xmin, ymax = -ymin;

	double temperature 					= 0.0;
	double baryon_chemical_potential 	= 0.0;
	double strange_chemical_potential 	= 0.0;
	double electric_chemical_potential 	= 0.0;
	double energy_density 				= 0.0;
	double baryon_density 				= 0.0;
	double strange_density 				= 0.0;
	double electric_density 			= 0.0;
	double entropy_density 				= 0.0;

	Vector<int,D> r0;
	r0.x[0] = 0.0;
	r0.x[1] = 0.0;

	string physical_quantities_filename = "./outputfiles/physical_quantities_timestep"
											+ to_string(countEP) + ".dat";
	ofstream physical_quantities_output(physical_quantities_filename.c_str());

	for (double x_local = xmin; x_local <= xmax + 1e-10; x_local += dx )
	for (double y_local = ymin; y_local <= ymax + 1e-10; y_local += dy )
	{
		r0.x[0] = x_local;
		r0.x[1] = y_local;
		compute_physical_quantities( linklist, r0, temperature,
			baryon_chemical_potential, strange_chemical_potential,
			electric_chemical_potential, energy_density, baryon_density,
			strange_density, electric_density, entropy_density );

		physical_quantities_output << setw(12) << setprecision(8) << scientific
			<< linklist.t << "   " << x_local << "   " << y_local << "   "
			<< temperature << "   " << baryon_chemical_potential << "   "
			<< strange_chemical_potential << "   " << electric_chemical_potential << "   "
			<< energy_density << "   " << baryon_density << "   " << strange_density << "   "
			<< electric_density << "   " << entropy_density << endl;
	}

	physical_quantities_output.close();

	return;
}




template <int D>
void Output<D>::compute_physical_quantities( LinkList<D> &linklist,
		Vector<int,D> r0, double & temperature, double & baryon_chemical_potential,
		double & strange_chemical_potential, double & electric_chemical_potential,
		double & energy_density, double & baryon_density, double & strange_density,
		double & electric_density, double & entropy_density )
{
	// define physical quantities to output to grid
	temperature 				= 0.0;
	baryon_chemical_potential 	= 0.0;
	strange_chemical_potential 	= 0.0;
	electric_chemical_potential = 0.0;
	energy_density 				= 0.0;
	baryon_density 				= 0.0;
	strange_density 			= 0.0;
	electric_density 			= 0.0;
	entropy_density 			= 0.0;

	double normalization 		= 0.0;

	// loop over SPH particles
	for (int iSPH = 0; iSPH < linklist.n(); iSPH++)
	{
		double kern 				 = kernel(r0-linklist._p[iSPH].r);
		normalization 				+= kern;
		energy_density 				+= kern * linklist._p[iSPH].EOSe();
		baryon_density 				+= kern * linklist._p[iSPH].EOSB();
		strange_density 			+= kern * linklist._p[iSPH].EOSS();
		electric_density 			+= kern * linklist._p[iSPH].EOSQ();
		temperature 				+= kern * linklist._p[iSPH].EOST();
		baryon_chemical_potential 	+= kern * linklist._p[iSPH].EOSmuB();
		strange_chemical_potential 	+= kern * linklist._p[iSPH].EOSmuS();
		electric_chemical_potential += kern * linklist._p[iSPH].EOSmuQ();
		entropy_density 			+= kern * linklist._p[iSPH].EOSs();
	}

	// normalize results
	energy_density 				/= normalization;
	baryon_density 				/= normalization;
	strange_density 			/= normalization;
	electric_density 			/= normalization;
	temperature 				/= normalization;
	baryon_chemical_potential 	/= normalization;
	strange_chemical_potential 	/= normalization;
	electric_chemical_potential /= normalization;
	entropy_density 			/= normalization;

	return;
}
//*/


template <int D>
void Output<D>::bsqsveprofile(LinkList<D> &linklist)
{
    ofstream EPN;

    countEP+=1;
    string epin ("bsqsveprofile");
    string dat (".dat");
    string sCEP;
    sCEP=convertInt(countEP);
    string epname;

    if (linklist.average==1)
    {
        epname=cfolder+"average_eprof"+ sCEP + ".dat";
    }
    else if (linklist.fcount==0&&linklist.average!=1)
    {
        epname=ofolder+epin+sCEP+dat;
    }
    else
    {
        string under ("_ev");
        string even;
        even=convertInt(linklist.fnum);
        epname=cfolder+epin+sCEP+under+even+dat;
    }

    //cout << epname << endl;
    EPN.open(epname.c_str());
    if (!EPN.is_open())
    {

        cout << "Error: cannot open eprofile" << sCEP <<".dat file!" << endl;
        exit(1);
    }
    else
    {
        EPN << linklist.t << endl;

        for (int i=0; i<linklist.n(); i++)
        {

            EPN << i << " " << linklist.t << " "
				<< linklist._p[i].r   << " "
				<< linklist._p[i].EOSp() << " "
				<< linklist._p[i].EOST()*197.3 << " "
				<< linklist._p[i].EOSmuB()*197.3 << " "
				<< linklist._p[i].EOSmuS()*197.3 << " "
				<< linklist._p[i].EOSmuQ()*197.3 << " "
				<< linklist._p[i].EOSe()*197.3 << " "
				<< linklist._p[i].EOSB() << " "
				<< linklist._p[i].EOSS() << " "
				<< linklist._p[i].EOSQ() << " "
				<< linklist._p[i].EOSs() << " "	//column 14
				<< linklist._p[i].eta/(linklist._p[i].gamma*linklist.t) << " "
				<< linklist._p[i].eta_sigma << " "
				<< linklist._p[i].sigma << " " 
				<< linklist._p[i].sigmaweight << " "
				<< linklist._p[i].stauRelax << " " 
				<< linklist._p[i].bigtheta   << " "
				<<  sqrt( linklist._p[i].shv.x[0][0]*linklist._p[i].shv.x[0][0]
						-2*linklist._p[i].shv.x[0][1]*linklist._p[i].shv.x[0][1]
						-2*linklist._p[i].shv.x[0][2]*linklist._p[i].shv.x[0][2]
						+ linklist._p[i].shv.x[1][1]*linklist._p[i].shv.x[1][1]
						+ linklist._p[i].shv.x[2][2]*linklist._p[i].shv.x[2][2]
						+2*linklist._p[i].shv.x[1][2]*linklist._p[i].shv.x[1][2]
						+pow(linklist.t,4)*linklist._p[i].shv33*linklist._p[i].shv33 ) << " "
				<< linklist._p[i].stauRelax/linklist.t * linklist._p[i].bigtheta << " "
				<< linklist._p[i].u.x[0]/linklist._p[i].gamma << " "
				<< linklist._p[i].u.x[1]/linklist._p[i].gamma << " "
				<< linklist._p[i].gamma << endl;

        }
    }

    EPN.close();

}
template <int D>
void Output<D>::eccout(LinkList<D> &linklist)
{
    ofstream EPN;

    countEP+=1;
    string epin ("eccCM");
    string dat (".dat");
    string sCEP;

    sCEP=convertInt(linklist.fnum);
    string epname;

    if (linklist.fcount==0&&linklist.average!=1)
    {
        epname=ofolder+epin+sCEP+dat;
    }
    else
    {
        string even=convertInt(linklist.fnum);
        epname=cfolder+epin+sCEP+dat;
    }


    EPN.open(epname.c_str());
    if (!EPN.is_open())
    {

        cout << "Error: cannot open eprofile" << sCEP <<".dat file!" << endl;
        exit(1);
    }

    double psi,rout;
    double e22=ecc(linklist,psi,rout,2,2);
    EPN << convertInt(linklist.fnum) << " " <<  e22 << " " << psi << " "  ;
    double e33=ecc(linklist,psi,rout,3,3);
    EPN << e33 << " " << psi  << " "  ;
    double e44=ecc(linklist,psi,rout,4,4);
    EPN << e44 << " " << psi  << " "  ;
    double e55=ecc(linklist,psi,rout,5,5);
    EPN << e55 << " " << psi  << " "  ;
    double e66=ecc(linklist,psi,rout,6,6);
    EPN << e66 << " " << psi  << " "  ;

    EPN << rout << endl;
    EPN.close();

}

template <int D>
double Output<D>::ecc(LinkList<D> &linklist, double  & psi,double & rout, int m,int n) {



    int max=linklist.n();

    double xcm=0,ycm=0,etot=0;
    for (int i=0; i<max; i++) {

        xcm+=linklist._p[i].r.x[0]*linklist._p[i].EOSe();
        ycm+=linklist._p[i].r.x[1]*linklist._p[i].EOSe();
        etot+=linklist._p[i].EOSe();

    }
    xcm/=etot;
    ycm/=etot;

    vector <double> r2,phi;
    r2.resize(max);
    phi.resize(max);
    double psit=0,psib=0,rb=0;
    for (int s=0; s<max; s++) {

        double xsub=(linklist._p[s].r.x[0]-xcm);
        double ysub=(linklist._p[s].r.x[1]-ycm);
        r2[s]=xsub*xsub+ysub*ysub;
        phi[s]=atan2(ysub,xsub);
        double rv=linklist._p[s].EOSe()*pow(r2[s],(m/2.));
        psit+=rv*sin(1.0*n*phi[s]);
        psib+=rv*cos(1.0*n*phi[s]);

        rb+=rv;

    }
    psit/=max;
    psib/=max;

    psi=1./(1.0*n)*atan2(psit,psib);
    double ec=0;
    for (int s=0; s<max; s++) ec+=linklist._p[s].EOSe()*pow(r2[s],m/2.)*cos(n*(phi[s]-psi));

    ec/=rb;
    rout=rb/etot;

    return ec;

}







template <int D>
string Output<D>::convertInt(int number)
{
    stringstream ss;//create a stringstream
    ss << number;//add number to the stream
    return ss.str();//return a string with the contents of the stream
}

template <int D>
void Output<D>::conservation(LinkList<D> &linklist)
{


    ESC.open(ESCname.c_str(), ios::out | ios::app );
//    ESC << linklist.t << " " <<   linklist.avgT*197.3 << " " << linklist.avgvt << " " <<  linklist.van <<   endl;
    ESC << linklist.t << " " <<  linklist.Eloss << endl;
    ESC.close();

    if (ESC.is_open())
    {
        cout << "error: still open" << endl;
        exit(1);
    }


}

//template <int D>
//void Output<D>::sconservation(LinkList<D> &linklist)
//{


//      ESC.open(ESCname.c_str(), ios::out | ios::app );
//    ESC << linklist.t << " " <<  linklist.Eloss << endl;
//      ESC.close();
//
//      if (ESC.is_open())
//      {
//          cout << "error: still open" << endl;
//         exit(1);
//      }
//
//
//}

template <int D>
void Output<D>::FOprint(LinkList<D> &linklist)
{

    FO.open(FOname.c_str(), ios::out | ios::app );
    for (int i=0; i< linklist.cf; i++)
    {
        FO <<  linklist.divTtemp[i] << " " << linklist.divT[i] << " " << linklist.gsub[i] << " " << linklist.uout[i] << " " << linklist.swsub[i] << " " <<   linklist.tlist[i] <<  " " << linklist.rsub[i] << " " << linklist.sFO[i] <<  " " << linklist.Tfluc[i] <<  endl;
    }
    FO.close();

    delete [] linklist.divTtemp;
    delete [] linklist.divT;
    delete [] linklist.gsub;
    delete [] linklist.uout;
    delete [] linklist.swsub;
    delete [] linklist.rsub;
    delete [] linklist.tlist;

    if (FO.is_open())
    {
        cout << "error: still open" << endl;
        exit(1);
    }


}

template <int D>
void Output<D>::vFOprint(LinkList<D> &linklist)
{

    FO.open(FOname.c_str(), ios::out | ios::app );
    for (int i=0; i< linklist.cf; i++)
    {
        FO <<  linklist.divTtemp[i] << " " << linklist.divT[i] << " " << linklist.gsub[i] << " " << linklist.uout[i] << " " << linklist.swsub[i] <<   " " << linklist.bulksub[i]  << " " <<   linklist.tlist[i] <<  " " << linklist.rsub[i] <<  " " << linklist.sFO[i] <<  " " << linklist.Tfluc[i] << endl;
    }
    FO.close();

    delete [] linklist.divTtemp;
    delete [] linklist.divT;
    delete [] linklist.gsub;
    delete [] linklist.uout;
    delete [] linklist.rsub;
    delete [] linklist.tlist;
    delete [] linklist.bulksub;

    if (FO.is_open())
    {
        cout << "error: still open" << endl;
        exit(1);
    }


}


template <int D>
void Output<D>::svFOprint(LinkList<D> &linklist)
{

    FO.open(FOname.c_str(), ios::out | ios::app );
    for (int i=0; i< linklist.cf; i++)
    {
        FO <<  linklist.divTtemp[i] << " " << linklist.divT[i] << " " << linklist.gsub[i] << " " << linklist.uout[i] << " " << linklist.swsub[i] << " " <<  linklist.bulksub[i] << " "  <<  linklist.shearsub[i].x[0][0] << " "  <<  linklist.shearsub[i].x[1][1] << " "  <<  linklist.shearsub[i].x[2][2] <<" " <<linklist.shear33sub[i] << " "  <<  linklist.shearsub[i].x[1][2] <<  " "  <<  linklist.tlist[i] <<  " " << linklist.rsub[i] <<  " " << linklist.sFO[i] <<   " " << linklist.Tfluc[i] <<endl; // removed entropy, added in tau
    }
    FO.close();

    delete [] linklist.divTtemp;
    delete [] linklist.divT;
    delete [] linklist.gsub;
    delete [] linklist.uout;
    delete [] linklist.swsub;
    delete [] linklist.bulksub;
    delete [] linklist.shearsub;
    delete [] linklist.rsub;
    delete [] linklist.tlist;
    delete [] linklist.shear33sub;

    if (FO.is_open())
    {
        cout << "error: still open" << endl;
        exit(1);
    }


}

template <int D>
void Output<D>::bsqsvFOprint(LinkList<D> &linklist)
{

    FO.open(FOname.c_str(), ios::out | ios::app );
    for (int i=0; i< linklist.cf; i++)
    {
        FO <<  linklist.divTtemp[i] << " " << linklist.divT[i] << " " << linklist.gsub[i] << " " << linklist.uout[i] << " " << linklist.swsub[i] << " " <<  linklist.bulksub[i] << " "  <<  linklist.shearsub[i].x[0][0] << " "  <<  linklist.shearsub[i].x[1][1] << " "  <<  linklist.shearsub[i].x[2][2] <<" " <<linklist.shear33sub[i] << " "  <<  linklist.shearsub[i].x[1][2] <<  " "  <<  linklist.tlist[i] <<  " " << linklist.rsub[i] <<  " " << linklist.sFO[i] <<   " " << linklist.Tfluc[i] <<endl; // removed entropy, added in tau
    }
    FO.close();

    delete [] linklist.divTtemp;
    delete [] linklist.divT;
    delete [] linklist.gsub;
    delete [] linklist.uout;
    delete [] linklist.swsub;
    delete [] linklist.bulksub;
    delete [] linklist.shearsub;
    delete [] linklist.rsub;
    delete [] linklist.tlist;
    delete [] linklist.shear33sub;

    if (FO.is_open())
    {
        cout << "error: still open" << endl;
        exit(1);
    }


}


template <int D>
void Output<D>::averages(LinkList<D> &linklist)
{

    ofstream AVG;
    string avgname=cfolder+"anisotropy.dat";

    AVG.open(avgname.c_str());
    if (!AVG.is_open())
    {

        cout << "Error: cannot open anisotropy.dat file!" << endl;
        exit(1);
    }
    else
    {
        for (int i=0; i<linklist.steps; i++)
        {
            AVG << linklist.tsave[i] << " " <<  linklist.pan[i] <<  " " <<  linklist.san[i] << endl;
        }
    }

    AVG.close();

    delete [] linklist.pan;
    delete [] linklist.san;

}

template <int D>
void Output<D>::cleardir(LinkList<D> &linklist)
{



    DIR *pdir;
    struct dirent *pent;
    if (linklist.fcount==0)
    {
        pdir=opendir(ofolder.c_str());
    }
    else
    {
        pdir=opendir(cfolder.c_str());
    }
    if (!pdir)
    {
        cout<<"Error: " << cfolder << " does not exist";
    }
    int errno3=0; //errno.h
    while ((pent=readdir(pdir)))
    {
        //cout<<pent->d_name;
        string file_delete=pent->d_name;
        file_delete= ofolder+file_delete;
//    DeleteFile(file_delete.c_str());
        std::remove(file_delete.c_str());
    }
    if (errno3)
    {
        cout<<"Error while accessing directory";
    }
    closedir(pdir);

}

template <int D>
void Output<D>::do_mkdir(const char *path, int mode)
{
    Stat            st;
//    int             status = 0;

//    cout << path << endl;
    if (stat(path, &st) != 0)
    {
        /* Directory does not exist */
        //  if (mkdir(path, mode) != 0)
        //    status = -1;
        cout << "Directory path doesn't exist" << endl;
    }
    else if (!S_ISDIR(st.st_mode))
    {
        errno = ENOTDIR;
        //  status = -1;
        cout << "Couldn't make the directory for output files" << endl;
        exit(1);
    }

}

#endif
