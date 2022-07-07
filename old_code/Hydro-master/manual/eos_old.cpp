#ifndef _EOS_CPP_
#define _EOS_CPP_

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
#include <string>
#include <stdio.h>
#include "eostables.h"

using namespace std;

//#include "fadiff.h"
//using namespace fadbad;


# include "newspline.hpp"
#include "tables.h"
#include "eos.h"
#include <float.h>

eos::eos()
{
	gg=3.297572;
	c=(16. +  7./2. *7.5 )/90. *cPI*cPI;
	cs=4./3.*pow(3.*c,0.25);
	cp=1./pow((4.*pow(c,0.25)),(4./3.));
	cp2=16./3.*cp;
	tempcutoff=65/197.3;
	ds1=-gg*pow(mpi,4)/(4.*cPI*cPI);
	es2=gg*pow(mpi,4)/(2.*cPI*cPI);
	es1=es2*3;
	ps1=gg/(2.*cPI*cPI);
	fourthirds=4./3.;
	onethird=1./3.;
	BIG=pow(1.,30.);
	
};

double eos::efreeze()
{
	temp=freezeoutT;
	double echeckout;
	if (typerem==table)
	{
		echeckout=spline(r8_T,y_e,temp,pos);
	}
	else if (typerem==ideal)
	{
		echeckout=3*c*pow(temp,4);
	}
	return echeckout;
	
};


double eos::sfreeze()
{
  temp=freezeoutT;
  double scheckout;
  if (typerem==table)
    {
      scheckout=spline(r8_T,y_s,temp,pos);
    }
  else if (typerem==ideal)
    {
      cout << "must setup still in EOS" << endl;
      exit(1);
    }
  return scheckout;
  
};

double eos::s_terms_T(double Tt)
{
	
	double etout,ptout;
	if (typerem==table)
	{
		if (Tt>ETH[0].T)
		{
		etout=spline(r8_T,y_e,Tt,pos);
		ptout=splinepos(y_T,y_p,Tt,pos);
		}
		else
		{
		etout=spline(r8_Tlow,y_elow,Tt,pos);
		ptout=splinepos(y_Tlow,y_plow,Tt,pos);
		}
	}
	else if (typerem==ideal)
	{
			
		ptout=c*pow(Tt,4);
		etout=3*ptout;
	}
	return (etout+ptout)/Tt;
	
};

double eos::cs2out(double Tt)
{
	double out;
	if (typerem==table)
	{
		double Ttsub=Tt+0.0001;
		double ee1,ee2,pp1,pp2;
		
		if (Tt<ETH[0].T)
		{
		ee1=spline(r8_Tlow,y_elow,Tt,pos);
		pp1=splinepos(y_Tlow,y_plow,Tt,pos);
		
		ee2=spline(r8_Tlow,y_elow,Ttsub,pos);
		pp2=splinepos(y_Tlow,y_plow,Ttsub,pos);
		
		}
		else 
		{
		ee1=spline(r8_T,y_e,Tt,pos);
		pp1=splinepos(y_T,y_p,Tt,pos);
		
		ee2=spline(r8_T,y_e,Ttsub,pos);
		pp2=splinepos(y_T,y_p,Ttsub,pos);
		
		}
		out=(pp2-pp1)/(ee2-ee1);
	}
	else if (typerem==ideal)
	{
		out=1./3.;
	}
	return out;
	
};

double eos::wfz(double Tt)
{
	double out;
	if (typerem==table)
	{
		out=spline(r8_T,y_e,Tt,pos)+spline(r8_T,y_p,Tt,pos);
	}
	else if (typerem==ideal)
	{
		out=4*c*pow(Tt,4.);
	}
	return out;
	
};


double eos::s_out(double eIC)
{
	double entro;
	if (typerem==ideal)
	{
		entro=cs*pow(eIC,0.75); 
	}
	else if (typerem==table)
	{
		
		if (eIC>echeck)
		{
			entro=spline(r8_e,y_s,eIC,pos);
		}
		else 
		{
			
			entro=spline(r8_elow,y_slow,eIC,pos);
		}
	}
		
	return entro;
}

double eos::e_out(double sIC)
{
  double ent;
  if (typerem==ideal)
    {
      ent=pow(sIC/cs,4./3.); 
    }
  else if (typerem==table)
    {
      
      if (sIC>scheck)
	{
	  ent=spline(r8_s,y_e,sIC,pos);
	}
      else 
	{
	  
	  ent=spline(r8_slow,y_elow,sIC,pos);
	}
    }
  
  return ent;
}

void eos::pressure_ideal()
{
	 pressure=pow( entropy,fourthirds)*cp;	 
	 
}


void eos::energy_ideal()
{

	energy=3*pressure;
};

void eos::temp_ideal()
{

	 temp=(energy+pressure)/entropy;
};

void eos::eosin(string type)
{
	
	typerem=type;
	
	ideal="ideal";
	table="table";
	
	
	
	
	if (type==table)
	{
		typerem=table;
		scheck=ETH[0].s+(ETH[1].s-ETH[0].s)/10000000;
		echeck=ETH[0].e+(ETH[1].e-ETH[0].e)/10000000;
	}
	else typerem="ideal";
	
};




void eos::temp_tab()
{

   temp=spline(r8_s,y_T,entropy,pos);
   
}

void eos::pressure_tab()
{
   pressure =splinepos(y_T,y_p,temp,pos);	
}

void eos::energy_tab()
{
   energy = splinepos(y_T,y_e,temp,pos);
};

double eos::Atable()
{
	Aout=w()-entropy*dwds();

	return Aout;
};

double eos::dwds_tab()
{
	
	double dwdstab;
	if (entropy<scheck) 
	{
	dwdstab=splinepos(y_Tlow,y_dtdslow,temp,pos);
	}
	else
	{ 
	dwdstab=splinepos(y_T,y_dtds,temp,pos);
	}
	
	
	
	//dtds is actually T+s*dtds
	return dwdstab;
};

void eos::outpos()
{
	
	cout << pos << " " << entropy << " " << scheck << endl;
	
};

double eos::dwds_ideal()
{
	return 4/3.*temp;
};


double eos::dwds()
{
	double dw_ds;
	
	if (typerem==ideal)
		dw_ds=dwds_ideal();
	else if (typerem==table)
	{
		dw_ds=dwds_tab();
	}
		
	return dw_ds;
	
};




eos::~eos()
{
	
	
};

double eos::e()
{
	return energy;
};

double eos::p()
{
	return pressure;
};

double eos::s()
{
	return entropy;
};

double eos::T()
{
	return temp;
};

double eos::A()
{

	if (typerem==ideal)
		return Aideal();
	else if (typerem==table)
		return Atable();
	
		
	return 0.;
};

double eos::Aideal()
{
	Aout=energy+pressure-cp2*pow(entropy,fourthirds);

	return Aout;
};



void eos::update_s(double s_in)
{
	entropy=s_in;

	
	if (typerem==ideal)
	{
		
		pressure_ideal();	
		energy_ideal();
		temp_ideal();
		
	}
	else if (typerem==table)
	{
		
		if (s_in<scheck)
		{
			
			tablow();
		}
		else
		{
		
		
		temp_tab();
		pressure_tab();
		energy_tab();
		}
		
	}
	epp();
	
};

void eos::tablow()
{
	temp=spline(r8_slow,y_Tlow,entropy,pos);

	pressure =splinepos(y_Tlow,y_plow,temp,pos);	
  	energy = splinepos(y_Tlow,y_elow,temp,pos);

};


double eos::dervs(double (eos::*f)(double) ,double x, double h,double & err)
{
 double out, con=1.4, con2=pow(con,2.), safe=2.;
 double errt, fac,hh, a[NTAB][NTAB];
 
 if (h==0.)
 {
 	cout << "Error: In double dervs - h must be nonzero" << endl;
 	exit(1);
 }
 hh=h; 
 a[0][0]=((this->*f)(x+hh)-(this->*f)(x-hh))/(2.*hh);
 err=BIG;
 
 for(int i=1;i<NTAB;i++)
 {
 	hh=hh/con;
 	a[0][i]=((this->*f)(x+hh)-(this->*f)(x-hh))/(2.*hh);
 	fac=con2;
 	
 	for (int j=1;j<=i;j++)
 	{
 		a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.);
 		fac=con2*fac;
 		errt=max(abs(a[j][i]-a[j-1][i]),abs(a[j][i]-a[j-1][i-1]));
 		if (errt<=err)
 		{
 			err=errt;
 			out=a[j][i];
 		}
 	}
 	if (abs(a[i][i]-a[i-1][i-1])>=(safe*err)) return out;
 }
 return out;

};

double eos::dervsprint(double (eos::*f)(double) ,double x, double h,double & err)
{
 
 double out, con=1.4, con2=pow(con,2.), safe=2.;
 double errt, fac,hh, a[NTAB][NTAB];
 
 if (h==0.)
 {
 	cout << "Error: In double dervs - h must be nonzero" << endl;
 	exit(1);
 }
 hh=h/2; 
 a[0][0]=((this->*f)(x+hh)-(this->*f)(x-hh))/(2.*hh);
 err=BIG;
 
 for(int i=1;i<NTAB;i++)
 {
 	hh=hh/con;
 	a[0][i]=((this->*f)(x+hh)-(this->*f)(x-hh))/(2.*hh);
 	fac=con2;
 	
 	for (int j=1;j<=i;j++)
 	{
 		a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.);
 		fac=con2*fac;
 		errt=max(abs(a[j][i]-a[j-1][i]),abs(a[j][i]-a[j-1][i-1]));
 		cout << i << " " << j << " " << a[0][i] << " " << a[j][i] << " " << err << endl;
 		if (errt<=err)
 		{
 			err=errt;
 			out=a[j][i];
 		}
 	}
 	if (abs(a[i][i]-a[i-1][i-1])>=(safe*err)) return out;
 }
 return out;

};


double eos::h_calc(double x)
{
	double ef;
	ef=DBL_EPSILON;
	
	if (x==0.) return 0.0001*pow(DBL_EPSILON,onethird);

 	return pow(ef,onethird)*x;

};

double eos::e_plow(double s)
{
	double e_sub,p_sub;
	
	
	e_sub = spline(r8_slow,y_elow,s,pos);
	p_sub = splinepos(y_slow,y_plow,s,pos);
	
	
	return e_sub+p_sub;
};

double eos::e_ptable(double s)
{
	double e_sub,p_sub;
	
	
	e_sub = spline(r8_s,y_e,s,pos);
	p_sub = splinepos(y_s,y_p,s,pos);
	
	
	return e_sub+p_sub;
};

#endif

