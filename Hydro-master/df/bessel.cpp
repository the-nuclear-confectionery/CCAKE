#ifndef _Bessel_CPP_
#define _Bessel_CPP_


#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <sstream>
#include <cmath>

#include "bessel.h"

double Bessel::K0(double x)
{

	double bes;
	if (x<2)
	{
		double y=x*x/4;
		bes=-log(x/2)*I0(x)+-0.57721566+y*(0.42278420+ y*(0.23069756+ y*(0.03488590+y*(0.00262698+ y*(0.00010750+y*0.0000074)))));
	}
	else
	{
		double y=2/x;
		bes=exp(-x)/sqrt(x)*(1.25331414+y*(-0.07832358+ y*(0.02189568+ y*(-0.01062446+y*(0.00587872+ y*(-0.00251540+y*0.00053208))))));
	}

	return bes;


}






double Bessel::I0(double x)
{

	double bes;
	double xp=abs(x);
	if (xp<3.75)
	{
		double y=pow(x/3.75,2);
		bes=1+y*(3.5156229+ y*(3.0899424+ y*(1.2067492+y*(0.2659732+ y*(0.0360768+y*0.0045813)))));
	}
	else
	{
		double y=3.75/xp;
		bes=exp(xp)/sqrt(xp)*(0.39894228+y*(0.01328592+ y*(0.00225319+ y*(-0.00157565+y*(0.00916281+ y*(-0.02057706+y*(0.02635537+y*(-0.01647633+0.00392377*y  )   )  ))))));
	}

	return bes;


}



      
double Bessel::K1(double x)
{

	double bes;
	double x1=x/2;
	if (x<2)
	{
		double y=x1*x1;
		bes=log(x1)*I1(x)+(1+y*(0.15443144+ y*(-0.67278579+ y*(-0.18156897+y*(-0.01919402+ y*(-0.00110404+y*-0.00004686))))))/x;
	}
	else
	{
		double y=1/x1;
		bes=exp(-x)/sqrt(x)*(1.25331414+y*(0.23498619+ y*(-0.03655620+ y*(0.01504268+y*(-0.00780353+ y*(0.00325614+y*-0.00068245))))));
	}

	return bes;


}





double Bessel::I1(double x)
{

	double bes;
	double xp=abs(x);
	if (xp<3.75)
	{
		double y=pow(x/3.75,2);
		bes=x*(0.5+ y*(0.87890594+ y*(0.51498869+y*(0.15084934+ y*(0.02658733+y*(0.00301532+ 0.00032411*y )  )))));
	}
	else
	{
		double y=3.75/xp;
		bes=exp(xp)/sqrt(xp)*(0.39894228+y*(-0.03988024+ y*(-0.00362018+ y*(0.00163801+y*(-0.01031555+ y*(0.02282967+y*(-0.02895312+y*(0.01787654-0.00420059*y  )   )  ))))));
		if (x<0) bes=-bes;
	}

	return bes;


}






double Bessel::Kn(int n, double x)
{

	
	if (n<2) {
	cout << "n too small in Bessel function" << endl;
	exit(1);
	}

	double tox=2/x;
	double bkm=K0(x);
	double bk=K1(x);
	for(int j=1;j<=(n-1);j++)
	{
		double bkp=bkm+j*tox*bk;
		bkm=bk;
		bk=bkp;

	}

	return bk;


}

      

#endif
