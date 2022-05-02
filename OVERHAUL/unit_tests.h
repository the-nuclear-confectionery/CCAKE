#include <cstdlib>
#include <iomanip>
#include <iostream>

//#include "matrix.h"
#include "vector.h"

#define RESET   "\033[0m"
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */

using namespace std;

inline void passed(const string & f) { cout << f << ": " << GREEN << "PASSED\n" << RESET; }
inline void failed(const string & f) { cout << f << ": " << RED << "FAILED\n" << RESET; }

template <class T, int D>
bool operator==( const Vector<T,D>& a, const Vector<T,D>& b )
{
  bool result = true;
  for (int i = 0; i < D; i++) result = result && ( a(i) == b(i) );
  return result;
}

////////////////////////////////////////////////////////////////////////////////
void check_copy()
{
  Vector<double, 2> a, b;
  a(0) = 1; a(1) = 2;
  b = a;
  if ( a == b ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
}

////////////////////////////////////////////////////////////////////////////////
void check_assignment()
{
  Vector<double, 2> a = 7;
  if ( a(0) == 7 && a(1) == 7 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
}



////////////////////////////////////////////////////////////////////////////////
void check_inplace_addition()
{
  Vector<double, 2> a = 7., b = 1.;
  a += b;
  if ( a(0) == 8 && a(1) == 8 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
}

////////////////////////////////////////////////////////////////////////////////
void check_inplace_subtraction()
{
  Vector<double, 2> a = 7., b = 1.;
  a -= b;
  if ( a(0) == 6 && a(1) == 6 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
}


////////////////////////////////////////////////////////////////////////////////
void check_inplace_multiplication()
{
  Vector<double, 2> a = 7.;
  double b = 3;
  a *= b;
  if ( a(0) == 21 && a(1) == 21 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
}

////////////////////////////////////////////////////////////////////////////////
void check_negation()
{
  Vector<double, 2> a = 7.;
  a = -a;
  if ( a(0) == -7 && a(1) == -7 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
}


////////////////////////////////////////////////////////////////////////////////
void check_addition()
{
  Vector<double, 2> a = 7, b = 3;
  Vector<double, 2> c = a + b;
  if ( c(0) == 10 && c(1) == 10 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
}


////////////////////////////////////////////////////////////////////////////////
void check_subtraction()
{
  Vector<double, 2> a = 7, b = 3;
  Vector<double, 2> c = a - b;
  if ( c(0) == 4 && c(1) == 4 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
}


////////////////////////////////////////////////////////////////////////////////
void check_multiplication()
{
  Vector<double, 2> a = 7;
  double b = 3;
  Vector<double, 2> c = b*a;
  if ( c(0) == 21 && c(1) == 21 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
}



////////////////////////////////////////////////////////////////////////////////
void check_inner()
{
  Vector<double, 2> a = 7;
  Vector<double, 2> b = 3;
  a(1) -= 1; b(0) += 1;
  if ( inner(a,b) == 46 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
}


////////////////////////////////////////////////////////////////////////////////
void check_Norm()
{
  Vector<double, 2> a = 3;
  a(1) += 1;
  if ( Norm(a) == 5 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
}




////////////////////////////////////////////////////////////////////////////////
void check_Norm2()
{
  Vector<double, 2> a = 3;
  a(1) += 1;
  if ( Norm2(a) == 25 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void run_vector_tests()
{
  check_copy();
  check_assignment();
  check_inplace_addition();
  check_inplace_subtraction();
  check_inplace_multiplication();
  check_negation();
  check_addition();
  check_subtraction();
  check_multiplication();
  check_inner();
  check_Norm();
  check_Norm2();
}