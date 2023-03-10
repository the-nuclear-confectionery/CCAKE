#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <vector>

#include "mathdef.h"
#ifdef DEBUG
#include <cassert>
#endif

using namespace std;


///\todo There are two things that I leave as suggestion for future work:
/// 1. We could replace the operator() with the [] operator, which is more
///    intuitive for the user. However, this would require to run thorough
///    the code and replace all the occurrences of operator() with [].
/// 2. There are quite a few instances of type casting in the class.
///    This is not a problem, but it would be nice to avoid it. It would also
///    make the class simpler by avoiding the use of template U and even 
///    possibly imporve build time.

///\todo We should also implement asserts on the dimension of the vectors. These
/// asserts should be turned off in release mode. A pratcical way to do this is
/// enclosing them in #ifdef DEBUG ... #endif
template <class T, int D>
class Vector
{
private:
//  vector<T> x{vector<T>(D,0)};
  T x[D];
public:
  Vector<T,D>(){for (int i = 0; i < D; i++) x[i] = (T)0.0;}
  Vector<T,D>(T x0){for (int i = 0; i < D; i++) x[i] = x0;}
  Vector<T,D>(std::initializer_list<T> x0){
    #ifdef DEBUG
    assert(x0.size() == D);
    #endif
    std::copy(x0.begin(), x0.end(), x);
  }


  template <class U> Vector<T,D>& operator=(const Vector<U,D>&);
  template <class U> Vector<T,D>& operator=(U);
  template <class U> Vector<T,D>& operator+=(const Vector<U,D>&);
  template <class U> Vector<T,D>& operator-=(const Vector<U,D>&);
  template <class U> Vector<T,D>& operator*=(U);
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

#include "../src/vector.cpp"

#endif