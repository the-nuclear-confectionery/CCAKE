#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <vector>

#include "mathdef.h"

using namespace std;

template <class T, int D>
class Vector
{
private:
//  vector<T> x{vector<T>(D,0)};
  T x[D];
public:
  Vector<T,D>(){}

  template <class U> Vector<T,D>& operator=(const Vector<U,D>&);
  Vector<T,D>& operator=(double);
  Vector<T,D>& operator+=(const Vector<T,D>&);
  Vector<T,D>& operator-=(const Vector<T,D>&);
  Vector<T,D>& operator*=(T);
  inline T& operator()(const int i)       { return x[i]; }
  inline T  operator()(const int i) const { return x[i]; }
};

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
template <class T, int D> double Norm2( const Vector<T,D>&              );

#include "vector.cpp"

#endif