#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <vector>

#include "SPH.h"
#include "int.h"
#include "spectra.h"
#include "tables.h"

using namespace std;



int main (int argc, char *argv[])
{
	_inputIC ics;
	
	if (argv[1])
	{
	ics.man= argv[1];
	
	
	
	if (argv[3])
	{
	stringstream s,s1;
	s << argv[3];
	s >> ics.start;
	
	ics.end=ics.start;
	ics.on=1;
	}
	else {	
	ics.on=0;
	}
	
	}
	else
	{
	ics.man="input.dat";
	ics.on=0;
	}

	





	SPH<2,3> sph;
	sph.readin(ics);
	int qmst=0;
	if ((sph.typ==1)||(sph.typ==3) )
  	{
  		sph.flist();	
  	}
  	
  	if (argv[2])
	{
	sph.rnum=argv[2];
	
	}
	if (argv[4])
	{
	sph.neg=argv[4];
	}
	if (argv[5])
	{
	stringstream s2;
	s2 << argv[5];
	s2 >> qmst;
	
	}

	

  	
  	
  	list l(sph.pt.size(),sph.phi.size());
  	l.setup(sph.pt,sph.phi);
  	
  	
  	

	
	string ofolder="out/"+sph.folder;
	string zero="0";
	if (argv[6]&&argv[6]!=zero)  {
           string newout=argv[6];
           ofolder=ofolder+"/"+newout;
        }
	if (qmst==1) ofolder=ofolder+"/qm";
	const double pisub=pow(2*PI,3);
	 //switches end result into 1/[GeV]^2
	
       

	if (sph.typ==1) //bulk
	{
		for (int i=0;i<sph.NHAD;i++) // sets up bulk coefficients for each hadron
		{
		sph.had[i].vfac=sph.had[i].deg/pisub;
		sph.setcoef(i);
		sph.calcsPI(i);
		sph.had[i].form3=0.25*pow(sph.had[i].mass,3);
		sph.had[i].halm2=0.5*pow(sph.had[i].mass,2);
		}
	}
	else if (sph.typ==2) //shear
	{
	for (int i=0;i<sph.NHAD;i++) // sets up bulk coefficients for each hadron
		{
		sph.had[i].form3=0.25*pow(sph.had[i].mass,3);
		sph.had[i].halm2=0.5*pow(sph.had[i].mass,2);
		sph.had[i].svfac=sph.had[i].deg/pisub*sph.scale;
		}
	}
	else if (sph.typ==3) //bulk+shear
	{
		for (int i=0;i<sph.NHAD;i++) // sets up bulk coefficients for each hadron
		{
		sph.had[i].vfac=sph.had[i].deg/pisub;
		sph.had[i].svfac=sph.had[i].vfac*sph.scale;
		sph.had[i].form3=0.25*pow(sph.had[i].mass,3);
		sph.had[i].halm2=0.5*pow(sph.had[i].mass,2);
		sph.setcoef(i);
		sph.calcsPI(i);
		}	
	}
	else {
		for (int i=0;i<sph.NHAD;i++) // sets up bulk coefficients for each hadron
		{
		sph.had[i].vfac=sph.had[i].deg/pisub;
		}
	}
	
	sph.printrun(ofolder);
  	cout << ofolder << endl;
	
	
	ofstream OUT;
	for (int ev=sph.start;ev<=sph.end;ev++) // runs over all the events
	{	
	
	  
	sph.readin2(ev);

	
	//multiply shear spectra prefactor by entropy
//	if (sph.typ>=2){
//	if (ev==sph.start){
//		for (int i=0;i<sph.NHAD;i++) 
//		{
//		sph.had[i].svfac/=sph.s;
//		}
//	}
//	}
	
	string conev=sph.convertInt(ev);
	string out1,out2;
	printstart<2,3>(ofolder,conev,out1,out2,sph);
	
	
	sph.checknu();
	
	
	
	for (int h=0;h<sph.NHAD;h++) // runs over the number of hadrons
	{
		
		
	  //		cout << "\r" << sph.had[h].id << " " << h << "= Hadron #"  <<flush;

	        if (sph.had[h].null!=2){

		if (sph.typ==0){
		for (int ps=0;ps<l.pTmax;ps++)
		{
			for(int i=0;i<l.phimax;i++) l.dNdpdphi.x[ps][i]=sph.dNdpdphi(l.pt.x[ps][i],l.phi.x[ps][i],sph.had[h]);
	 	}
	 	}
	 	else{
	 	for (int ps=0;ps<l.pTmax;ps++)
		{
			for(int i=0;i<l.phimax;i++) {
			l.dNdpdphi.x[ps][i]=sph.dNdpdphi(l.pt.x[ps][i],l.phi.x[ps][i],sph.had[h]);
			l.dNdpdphic.x[ps][i]=sph.outc;
//			if (l.dNdpdphi.x[ps][i]<0||isnan(l.dNdpdphi.x[ps][i])) l.dNdpdphi.x[ps][i]=0;
//			if (l.dNdpdphic.x[ps][i]<0||isnan(l.dNdpdphic.x[ps][i])) l.dNdpdphic.x[ps][i]=0;
			if (isnan(l.dNdpdphi.x[ps][i])) l.dNdpdphi.x[ps][i]=0;
			if (isnan(l.dNdpdphic.x[ps][i])) l.dNdpdphic.x[ps][i]=0;
			}
	 	}
	 	
	 	}
		
		
		//print off spectra
		if (sph.typ>0) printc(out2,h,l,sph);
//		else print(out1,h,l,sph);
 		print(out1,h,l,sph);
		
		if (sph.had[h].anti==1){
		  int nsub=-sph.had[h].id;
		  if (sph.typ>0) nprintc(out2,h,l,sph,nsub);
		  else nprint(out1,h,l,sph,nsub);
		  
		}

                if (sph.had[h].null==1){
		  int ssid=sph.had[h].sids.size();
		  for (int ss=0;ss<ssid;ss++){
		    int hsub=sph.had[h].sids[ss];
		    if (sph.typ>0) nprintc(out2,h,l,sph,sph.had[hsub].id);
		    else nprint(out1,h,l,sph,sph.had[hsub].id);
		    
		    if (sph.had[hsub].anti==1){
		      if (sph.typ>0) nprintc(out2,h,l,sph,-sph.had[hsub].id);
		      else nprint(out1,h,l,sph,-sph.had[hsub].id);
		    }
		  }
		}


		
		}
		
	 	
		
		
		
	}
	
	
	delete [] sph.par;
	delete [] sph.qv;
	}
	
	if (sph.typ>0){
	 l.destroyc();
	}
	else l.destroy();
	
	
	if (sph.typ>0) delete [] sph.hlist;
	 

 return 0;
}
