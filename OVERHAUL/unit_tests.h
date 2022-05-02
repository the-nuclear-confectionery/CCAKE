#include <cstdlib>
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

void check_assignment()
{
  Vector<int, 2> a, b;
  a(0) = 1; a(1) = 2;
  b = a;
  if ( a == b ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
}

void run_vector_tests()
{
  check_assignment();
}


/*
////////////////////////////////////////////////////////////////////////////////
void run_vector_and_matrix_unit_tests()
{

  template <class U> Vector<T,D>& operator=(const Vector<U,D>&);
  Vector<T,D>& operator=(double);
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