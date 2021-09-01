#include <bits/stdc++.h>
#include <cmath>
#include <iomanip>

#define EPSILON 0.0000001

using namespace std;

double func(double x)
{
	return cos(x) - x*x;
}

double derivFunc(double x)
{
	return -sin(x) - 2.0*x;
}

void newtonRaphson(double x)
{
	double h = func(x) / derivFunc(x);
	while (abs(h) >= EPSILON)
	{
		h = func(x)/derivFunc(x);

		// x(i+1) = x(i) - f(x) / f'(x)
		x = x - h;
	}

	cout << "The value of the root is : " << setprecision(12) << x << endl;
}

double xmap(double t, double a, double b)
{
	return a + (b-a)/(1.0+exp(-t));
}

double dxmap_dt(double t, double a, double b)
{
	double et = exp(t);
	return (b-a)*et/(1.0+et*et);
}

void newtonRaphson(double a, double b)
{
	double t = 0.0;
	double h = func(xmap(t,a,b)) / (derivFunc(xmap(t,a,b))*dxmap_dt(t,a,b));
	while (abs(h) >= EPSILON)
	{
		h = func(xmap(t,a,b)) / (derivFunc(xmap(t,a,b))*dxmap_dt(t,a,b));
		t -= h;
	}

	cout << "The value of the root is : t = " << setprecision(12) << t
		<< " or x = " << xmap(t,a,b) << endl;
}


// Driver program to test above
int main()
{
	double x0 = 20; // Initial values assumed
	newtonRaphson(x0);
	newtonRaphson(0.0, 2.0);
	return 0;
}

