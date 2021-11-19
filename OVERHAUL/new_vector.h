#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <vector>

#include "mathdef.h"
#include "vector.h"

using namespace std;

template <class T, int D>
class Vector
{
private:
  vector<T> x(D,0);
public:
  Vector<T,D>(){}

  template <class U> Vector<T,D>& operator=(Vector<U,D>);
  Vector<T,D>& operator=(double);
  Vector<T,D>& operator+=(Vector<T,D>);
  Vector<T,D>& operator-=(Vector<T,D>);
  Vector<T,D>& operator*=(T);
  T operator[](int i) { return x[i]; }
};

//==============================================================================
// overloaded operations
template <class T, int D> Vector<T,D> operator+ ( Vector<T,D>, Vector<T,D> );
template <class T, int D> Vector<T,D> operator- ( Vector<T,D>              );
template <class T, int D> Vector<T,D> operator- ( Vector<T,D>, Vector<T,D> );
template <class T, int D> Vector<T,D> operator* ( T, Vector<T,D>           );
template <class T, int D> ostream&    operator<<( ostream&, Vector<T,D>    );

//==============================================================================
// vector functions
template <class T, int D> double inner( Vector<T,D>, Vector<T,D> );
template <class T, int D> double norm(  Vector<T,D>              );
template <class T, int D> double norm2( Vector<T,D>              );

#include "vector.cpp"

#endif