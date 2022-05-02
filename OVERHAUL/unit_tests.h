//#include "matrix.h"
#include "vector.h"

using namespace std;

template <class T, int D>
bool operator==( const Vector<T,D>& a, const Vector<T,D>& b )
{
  bool result = true;
  for (int i = 0; i < D; i++) result = result && ( a(i) == b(i) );
  return result;
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

*/

}