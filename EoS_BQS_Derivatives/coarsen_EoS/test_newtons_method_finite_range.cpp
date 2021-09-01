#include <bits/stdc++.h>
#include <cmath>

#define EPSILON 0.001

using namespace std;

double func(double x)
{
	return cos(x) - x*x;
}

double derivFunc(double x)
{
	return -sin(x) - 2.0*x;
}

// Function to find the root
void newtonRaphson(double x)
{
	double h = func(x) / derivFunc(x);
	while (abs(h) >= EPSILON)
	{
		h = func(x)/derivFunc(x);

		// x(i+1) = x(i) - f(x) / f'(x)
		x = x - h;
	}

	cout << "The value of the root is : " << x << endl;
}

// Driver program to test above
int main()
{
	double x0 = 20; // Initial values assumed
	newtonRaphson(x0);
	return 0;
}

