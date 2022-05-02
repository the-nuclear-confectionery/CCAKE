#include <cstdlib>
#include <iomanip>
#include <iostream>

//#include "matrix.h"
#include "vector.h"

#define RESET   "\033[0m"
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */

using namespace std;

inline void passed(const string & f) { cout << f << ":" << right << setw(50) << GREEN << "PASSED\n" << RESET; }
inline void failed(const string & f) { cout << f << ":" << right << setw(50) << RED << "FAILED\n" << RESET; }

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




/*
////////////////////////////////////////////////////////////////////////////////
void run_vector_and_matrix_unit_tests()
{

  Vector<T,D>& operator+=(const Vector<T,D>&);
  Vector<T,D>& operator-=(const Vector<T,D>&);
  Vector<T,D>& operator*=(T);


//==============================================================================
// overloaded operations
template <class T, int D> Vector<T,D> operator+ ( const Vector<T,D>&, const Vector<T,D>& );
template <class T, int D> Vector<T,D> operator- ( const Vector<T,D>&              );
template <class T, int D> Vector<T,D> operator- ( const Vector<T,D>&, const Vector<T,D>& );
template <class T, int D> Vector<T,D> operator* ( T, const Vector<T,D>&           );
template <class T, int D> ostream&    operator<<( ostream&, const Vector<T,D>&    );

//==============================================================================
// vector functions
template <class T, int D> double inner( const Vector<T,D>&, const Vector<T,D>& );
template <class T, int D> double Norm(  const Vector<T,D>&              );
template <class T, int D> double Norm2( const Vector<T,D>&              );  }



Vector<T,D>& Vector<T,D>::operator=(const Vector<U,D> & a)
Vector<T,D>& Vector<T,D>::operator=(double a)
Vector<T,D>& Vector<T,D>::operator+=(const Vector<T,D> & a)
Vector<T,D>& Vector<T,D>::operator-=(const Vector<T,D> & a)
Vector<T,D>& Vector<T,D>::operator*=(T l)
Vector<T,D> operator+(const Vector<T,D>& a, const Vector<T,D>& b)
Vector<T,D> operator-(const Vector<T,D>& a)
Vector<T,D> operator-(const Vector<T,D>& a, const Vector<T,D>& b)
Vector<T,D> operator*(T l, const Vector<T,D>& a)
double Norm(const Vector<T,D>& a)
double Norm2(const Vector<T,D>& a)
ostream& operator<<(ostream& os, const Vector<T,D>& a)
double inner (const Vector<T,D>& a, const Vector<T,D>& b)


}


*/



void run_vector_tests()
{
  check_copy();
  check_assignment();
  check_inplace_addition();
  check_inplace_subtraction();
  check_inplace_multiplication();
  check_negation();
}