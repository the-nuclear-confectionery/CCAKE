#ifndef MATRIX_CPP
#define MATRIX_CPP

#include <cmath>
#include <cstdio>
#include <type_traits>

#include "../include/matrix.h"
#include "../include/vector.h"

// D1 is the number for the rows, D2 is the number for the columns
template <class T, int D1, int D2>
template <class U>
Matrix<T,D1,D2>& Matrix<T,D1,D2>::operator=(const Matrix<U,D1,D2>& a)
{
  int k = 0;
  for (int i=0; i<D1; i++)
  for (int j=0; j<D2; j++)
    x[k++] = (T)a(i,j);
  return *this; // AOK
}


template <class T, int D1, int D2>
Matrix<T,D1,D2>& Matrix<T,D1,D2>::operator=(double a)
{
  int k = 0;
  for (int i=0; i<D1; i++)
  for (int j=0; j<D2; j++)
    x[k++] = (T)a;
  return *this; // AOK
}


template <class T, int D1, int D2>
Matrix<T,D1,D2>& Matrix<T,D1,D2>::operator+=(const Matrix<T,D1,D2>& a)
{
  int k = 0;
  for(int i=0; i<D1; i++)
  for(int j=0; j<D2; j++)
    x[k++] += a(i,j);
  return *this; // AOK
}


template <class T, int D1, int D2>
Matrix<T,D1,D2>& Matrix<T,D1,D2>::operator-=(const Matrix<T,D1,D2>& a)
{
  int k = 0;
  for(int i=0; i<D1; i++)
  for(int j=0; j<D2; j++)
    x[k++] -= a(i,j);
  return *this; // AOK
}


template <class T, int D1, int D2>
Matrix<T,D1,D2>& Matrix<T,D1,D2>::operator*=(T l)
{
  for (int i=0; i<D1*D2; i++) x[i] *= l;
  return *this; // AOK
}

template <class T, int D1, int D2>
Matrix<T,D1,D2>& Matrix<T,D1,D2>::operator*=(const Matrix<T,D2,D1>& b)
{
  //only works for matrices with the same dimensions
  static_assert(D1 == D2);

  if (true)
  {
    cout << "You should not be using this function!  This is buggy!" << endl;
    exit(1);
  }

//  for (int i=0; i<D1; i++)
//  for (int j=0; j<D1; j++)
//  {
//    T sub = 0.0;
//    for (int k=0; k<D2; k++) sub += x[index(i,k)]*b(k,j);
//    x[index(i,j)] = sub;
//  }
  x = x*b;
  return *this;
} // AOK


template <class T, int D1, int D2>
Matrix<T,D1,D2> operator+(const Matrix<T,D1,D2>& a, const Matrix<T,D1,D2>& b)
{
  Matrix<T,D1,D2> t;
  for (int i=0; i<D1; i++)
  for (int j=0; j<D2; j++)
    t(i,j) = a(i,j)+b(i,j);
  return t;
} // AOK


template <class T, int D1, int D2>
Matrix<T,D1,D2> operator-(const Matrix<T,D1,D2>& a)
{
  return (-1.0)*a;
} // AOK


template <class T, int D1, int D2>
Matrix<T,D1,D2> operator-(const Matrix<T,D1,D2>& a, const Matrix<T,D1,D2>& b)
{
  Matrix<T,D1,D2> t;
  for (int i = 0; i < D1; i++)
  for (int j = 0; j < D2; j++)
    t(i,j) = a(i,j) - b(i,j);
  return t;
} // AOK


template <class T, int D1, int D2>
Matrix<T,D1,D2> operator*(T l, const Matrix<T,D1,D2>& a)
{
  Matrix<T,D1,D2> t;
  for (int i = 0; i < D1; i++)
  for (int j = 0; j < D2; j++)
    t(i,j) = l*a(i,j);

  return t;
} // AOK

template <class T, int D1, int D2>
Matrix<T,D1,D2>& Matrix<T,D1,D2>::identity()
{
  static_assert( D1 == D2 );
//  if (D1!=D2) cout << "Error: not true identity matrix!" << endl;

  for (int i=0; i<D1; i++)
  for (int j=0; j<D2; j++)
    x[index(i,j)] = ( i == j ) ? 1 : 0;

  return *this;
} // AOK



template <class T, int D1, int D2>
ostream& operator<<(ostream& os, const Matrix<T,D1,D2>& a)
{
  for (int i=0; i<D1; i++)
  for (int j=0; j<D2; j++)
    os << a(i,j) << " ";
  return os;
} // AOK

template <class T, int D1, int D2>
Matrix<T,D2,D1> transpose(const Matrix<T,D1,D2>& a)
{
  Matrix<T,D2,D1> t;

  for (int i=0; i<D2; i++)
  for (int j=0; j<D1; j++)
    t(i,j) = a(j,i);

  return t;
} // AOK

template <class T, int D1, int D2, int Da2, int Db1>
Matrix<T,D1,D2> operator*(const Matrix<T,D1,Da2>& a, const Matrix<T,Db1,D2>& b)
{
  static_assert( Da2 == Db1 );

  Matrix<T,D1,D2> t;
  for (int i=0; i<D1; i++)
  for (int j=0; j<D2; j++)
  for (int k=0; k<Da2; k++)
    t(i,j) += a(i,k) * b(k,j);

  return t;
} // AOK

template <class T, int D1, int D2>
Vector<T,D1> operator*(const Matrix<T,D1,D2>& a, const Vector<T,D2>& b)
{
  Vector<T,D1> t;
  for (int i=0; i<D1; i++)
  for (int j=0; j<D2; j++)
    t(i) += a(i,j) * b(j);

  return t;
} // AOK


//template <class T, int D1, int D2>
//Matrix<T,D1,D2> outer( const Vector<T,D1>& a, const Vector<T,D2>& b )
template <class T, int D1, int D2>
Matrix<T,D1,D2> operator* ( const Vector<T,D1>& a, const Vector<T,D2>& b )
{
  Matrix<T,D1,D2> t;
  for (int i = 0; i < D1; i++)
  for (int j = 0; j < D2; j++)
    t(i,j) = a(i)*b(j);
  return t;
} // AOK


// takes spatial components of l^{th} row in space-time tensor (Matrix a)
template <class T, int D1, int D2>
Vector<T,(D2-1)> rowp1(int l, const Matrix<T,D1,D2>& a)
{
  Vector<T,(D2-1)> v;
  for (int i=1; i<D2; i++) v(i-1)=(T)a(l,i);
  return v;
} // AOK

// takes spatial components of l^{th} column in space-time tensor (Matrix a)
template <class T, int D1, int D2>
Vector<T,(D1-1)> colp1(int l, const Matrix<T,D1,D2>& a)
{
  Vector<T,(D1-1)> v;
  for(int i=1; i<D1; i++) v(i-1)=(T)a(i,l);
  return v;
} // AOK

template <class T, int D1, int D2>
Vector<T,D2> operator*(const Vector<T,D1>& a, const Matrix<T,D1,D2>& b)
{
  Vector<T,D2> t;
  for (int i=0; i<D1; i++)
  for (int j=0; j<D2; j++)
    t(j)+=a(i)*b(i,j);

  return t;
} // AOK


template <class T, int D1>
double deter(const Matrix<T,D1,D1>& a)
{
  // these are the only supported options for right now
  static_assert( D1==2 || D1==3 );
  if (D1==2)
    return a(0,0)*a(1,1)-a(0,1)*a(1,0);
  else if (D1==3)
    return   a(0,0)*(a(1,1)*a(2,2)-a(1,2)*a(2,1))
           + a(0,1)*(a(1,2)*a(2,0)-a(1,0)*a(2,2))
           + a(0,2)*(a(1,0)*a(2,1)-a(1,1)*a(2,0));
} // AOK

template <class T, int D1, int D2>
double con( const Matrix<T,D1,D2>& a, const Matrix<T,D1,D2>& b )
{
  double t = 0.0;
  for (int i = 0; i < D2; i++)
  for (int j = 0; j < D1; j++)
    t += a(j,i)*b(j,i);
  return t;
} // AOK

template <class T, int D1, int D2>
double con2(const Matrix<T,D1,D2>& a, const Matrix<T,D1,D2>& b)
{
  double t = 0.0;
  for (int i = 0; i < D2; i++)
  for (int j = 0; j < D1; j++)
    t += a(i,j)*b(i,j);
  return t;
} // AOK


// store transverse spatial components of a in b
template <class T, int D1, int D2>
void mini( Matrix<T,D1-1,D2-1> &b, const Matrix<T,D1,D2>& a)
{
  for(int i=1; i<D1; i++)
  for(int j=1; j<D2; j++)
    b(i-1,j-1) = (T)a(i,j);
} // AOK


// store a in transverse spatial components of b
template <class T, int D1, int D2>
void tmini( Matrix<T,D1,D2>& b, const Matrix<T,D1-1,D2-1>& a )
{
  for (int j = 0; j < (D2-1); j++)
  for (int i = 0; i < (D1-1); i++)
    b(i+1,j+1) = (T)a(i,j);
} // AOK


#endif