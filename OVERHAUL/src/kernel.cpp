#include <iomanip>

#include "../include/constants.h"
#include "../include/kernel.h"

using namespace constants;

using std::cout;
using std::endl;
using std::setprecision;

namespace kernel
{

  double knorm, knorm2, kgrad, kgrad2;

  void set_kernel_parameters( double h )
  {
//    knorm       = 10/7./pi/(_h*_h);
//    kgrad       = -10.0/7./pi/pow(_h,3)*3/4.;
//    kgrad2      = 10.0/7./pi/pow(_h,3)/_h;
    knorm       = 10.0/(7.*pi*h*h);
    knorm2      = knorm*0.25;
    kgrad       = -15.0/(14.0*pi*h*h*h);
    kgrad2      = 10.0/(7.0*pi*h*h*h*h);
//cout << "CHECK KERNEL: " << setprecision(16) << knorm << "   " << knorm2 << "   " << kgrad << "   " << kgrad2 << "   " << pi << endl;
//if (1) exit(-1);

    return;
  }

  double kernel( double q )
  {
    if ( q >= 2.0 )
      return 0.0;
    if ( q >= 1.0 )
      return knorm2*(2.0-q)*(2.0-q)*(2.0-q);

    double qq=q*q;
    return knorm*(1 - 1.5*qq + 0.75*q*qq);
  }


  double kernel( const Vector<double,2> & a, double h )
  {
    double q = Norm(a)/h;

    if ( q >= 2.0 )
      return 0.0;
    if ( q >= 1.0 )
      return knorm2*(2.0-q)*(2.0-q)*(2.0-q);

    double qq=q*q;
    return knorm*(1 - 1.5*qq + 0.75*q*qq);
  }


  Vector<double,2> gradKernel( const Vector<double,2> & a, double r, double h )
  {
    double q = r/h;

    if ( q >= 2.0 )
      return Vector<double,2>();
    if ( q >= 1.0 )
      return (kgrad/r)*(2.0-q)*(2.0-q)*a;

    return kgrad2*( -3.0+2.25*q )*a;
  }


  Vector<double,2> gradKernel( const Vector<double,2> & a, double h )
  {
    double r = Norm(a);
    double q = r/h;

    if ( q >= 2.0 )
      return Vector<double,2>();
    if ( q >= 1.0 )
      return (kgrad/r)*(2.0-q)*(2.0-q)*a;

    return kgrad2*( -3.0+2.25*q )*a;
  }

}