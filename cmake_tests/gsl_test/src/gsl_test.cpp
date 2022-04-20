#include <cstdio>
#include <cmath>
#include <iostream>
#include <vector>

#include <gsl/gsl_sf_bessel.h>

double f(double x) {return 1.0/(1.0+std::exp(x));}

int main (void)
{
  double x = 5.0;
  double y = gsl_sf_bessel_J0(x);
  printf ("J0(%g) = %.18e\n", x, y);

  long n = 100000000;
  for (long i = 1; i <= n; i++)
    x = f(x);

  std::cout << x << std::endl;

  return 0;
}
