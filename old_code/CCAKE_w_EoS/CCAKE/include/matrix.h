#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>

#include "mathdef.h"
#include "vector.h"

using namespace std;

template <class T, int D1, int D2>
class Matrix
{
private:
//  vector<T> x{vector<T>(D1*D2,0)};
  T x[D1*D2];
  inline int index(const int i, const int j) const { return i*D2+j; }
public:
  Matrix<T,D1,D2>() {for (int i=0; i<D1*D2; i++) x[i] = 0;}
  Matrix<T,D1,D2>(T x0) {for (int i = 0; i < D1*D2; i++) x[i] = x0;}

  Matrix<T,D1,D2>& operator=(double);
  Matrix<T,D1,D2>& operator+=(const Matrix<T,D1,D2>&);
  Matrix<T,D1,D2>& operator-=(const Matrix<T,D1,D2>&);
  Matrix<T,D1,D2>& operator*=(const Matrix<T,D2,D1>& b);
  Matrix<T,D1,D2>& operator*=(T);
  Matrix<T,D1,D2>& identity();
  template <class U>
    Matrix<T,D1,D2>& operator=(const Matrix<U,D1,D2>& a);
  inline T& operator()(const int i, const int j)       { return x[index(i,j)]; }
  inline T  operator()(const int i, const int j) const { return x[index(i,j)]; }
};

//==============================================================================
// overloaded operator functions
//==============================
// print matrix
template <class T, int D1, int D2>
  ostream&        operator<<( ostream& os, const Matrix<T,D1,D2>& a );
//==============================
// sum two matrices
template <class T, int D1, int D2>
  Matrix<T,D1,D2> operator+ ( const Matrix<T,D1,D2>&, const Matrix<T,D1,D2>& );
//==============================
// negative of matrix
template <class T, int D1, int D2>
  Matrix<T,D1,D2> operator- ( const Matrix<T,D1,D2>& );
//==============================
// subtract two matrices
template <class T, int D1, int D2>
  Matrix<T,D1,D2> operator- ( const Matrix<T,D1,D2>&, const Matrix<T,D1,D2>& );
//==============================
// multiply matrix by scalar
template <class T, int D1, int D2>
  Matrix<T,D1,D2> operator* ( T, const Matrix<T,D1,D2>& );
//==============================
// standard matrix multiplication
template <class T, int D1, int D2, int Da2, int Db1>
  Matrix<T,D1,D2> operator* ( const Matrix<T,D1,Da2>& a, const Matrix<T,Db1,D2>& b );
//==============================
// left-multiply matrix a on vector b: a.b
template <class T, int D1, int D2>
  Vector<T,D1>    operator* ( const Matrix<T,D1,D2>& a, const Vector<T,D2>& b );
//==============================
// right-multiply matrix a on vector b: b.a
template <class T, int D1, int D2>
  Vector<T,D2>    operator* ( const Vector<T,D1>& a, const Matrix<T,D1,D2>& b );

//==============================================================================
// matrix functions
//==============================
template <class T, int D1, int D2>
  Matrix<T,D2,D1> transpose( const Matrix<T, D1, D2>& a );

//==============================
template <class T, int D1, int D2>
  Vector<T,(D2-1)> rowp1( int l, const Matrix<T, D1, D2>& a );
//==============================
template <class T, int D1, int D2>
  Vector<T,(D1-1)> colp1( int l, const Matrix<T, D1, D2>& a );

//==============================
template <class T, int D1>
  double deter( const Matrix<T,D1,D1>& a );
//==============================
template <class T, int D1, int D2>
  double con(   const Matrix<T,D1,D2>& a, const Matrix<T,D1,D2>& b );
//==============================
template <class T, int D1, int D2>
  double con2(  const Matrix<T,D1,D2>& a, const Matrix<T,D1,D2>& b );
//==============================
template <class T, int D1, int D2>
  void mini(    Matrix<T,D1-1,D2-1> &b, const Matrix<T,D1,D2>& a );
//==============================
template <class T, int D1, int D2>
  void tmini(   Matrix<T,D1,D2>& b, const Matrix<T,D1-1,D2-1>& a );

#include "../src/matrix.cpp"

#endif