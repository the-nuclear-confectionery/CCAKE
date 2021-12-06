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
  vector<vector<T> > x;
public:
  Matrix<T,D1,D2>();
  Matrix<T,D1,D2>& operator=(double);
  Matrix<T,D1,D2>& operator+=(Matrix<T,D1,D2>);
  Matrix<T,D1,D2>& operator-=(Matrix<T,D1,D2>);
  Matrix<T,D1,D2>& operator*=(Matrix<T,D2,D1> b);
  Matrix<T,D1,D2>& operator*=(T);
  Matrix<T,D1,D2>& identity();
  template <class U>
    Matrix<T,D1,D2>& operator=(Matrix<U,D1,D2> a);
  vector<T>& operator[](int i) { return x[i]; }
};

//==============================================================================
// overloaded operator functions
//==============================
// print matrix
template <class T, int D1, int D2>
  ostream&        operator<<( ostream& os, Matrix<T,D1,D2> a );
//==============================
// sum two matrices
template <class T, int D1, int D2>
  Matrix<T,D1,D2> operator+ ( Matrix<T,D1,D2>, Matrix<T,D1,D2> );
//==============================
// negative of matrix
template <class T, int D1, int D2>
  Matrix<T,D1,D2> operator- ( Matrix<T,D1,D2> );
//==============================
// subtract two matrices
template <class T, int D1, int D2>
  Matrix<T,D1,D2> operator- ( Matrix<T,D1,D2>, Matrix<T,D1,D2> );
//==============================
// multiply matrix by scalar
template <class T, int D1, int D2>
  Matrix<T,D1,D2> operator* ( T, Matrix<T,D1,D2> );
//==============================
// standard matrix multiplication
template <class T, int D1, int D2, int Da2, int Db1>
  Matrix<T,D1,D2> operator* ( Matrix<T,D1,Da2> a, Matrix<T,Db1,D2> b );
//==============================
// left-multiply matrix a on vector b: a.b
template <class T, int D1, int D2>
  Vector<T,D1>    operator* ( Matrix<T,D1,D2> a, Vector<T,D2> b );
//==============================
// right-multiply matrix a on vector b: b.a
template <class T, int D1, int D2>
  Vector<T,D2>    operator* ( Vector<T,D1> a, Matrix<T,D1,D2> b );

//==============================================================================
// matrix functions
//==============================
template <class T, int D1, int D2>
  Matrix<T,D2,D1> transpose( Matrix<T,D1,D1> a );

//==============================
template <class T, int D1, int D2>
  Vector<T,(D2-1)> rowp1( Matrix<T,D1,D1> a );
//==============================
template <class T, int D1, int D2>
  Vector<T,(D1-1)> colp1( Matrix<T,D1,D1> a );

//==============================
template <class T, int D1>
  double deter( Matrix<T,D1,D1> a );
//==============================
template <class T, int D1, int D2>
  double con(   Matrix<T,D1,D2> a, Matrix<T,D1,D2> b );
//==============================
template <class T, int D1, int D2>
  double con2(  Matrix<T,D1,D2> a, Matrix<T,D1,D2> b );
//==============================
template <class T, int D1, int D2>
  void mini(    Matrix<T,D1-1,D2-1> &b, Matrix<T,D1,D2> a );
//==============================
template <class T, int D1, int D2>
  void tmini(   Matrix<T,D1,D2> &b, Matrix<T,D1-1,D2-1>a );

#include "matrix.cpp"

#endif