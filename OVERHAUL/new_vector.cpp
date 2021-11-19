#ifndef VECTOR_CPP
#define VECTOR_CPP

#include <cmath>

#include "vector.h"

//==============================================================================
template <class T, int D>
template <class U>
Vector<T,D>& Vector<T,D>::operator=(Vector<U,D> a)
{
  for (int i = 0; i < D; i++) x[i] = (T)a[i];
  return *this;
}

//==============================================================================
template <class T, int D>
Vector<T,D>& Vector<T,D>::operator=(double a)
{
  for (int i = 0; i < D; i++) x[i] = (T)a;
  return *this;
}

//==============================================================================
template <class T, int D>
Vector<T,D>& Vector<T,D>::operator+=(Vector<T,D> a)
{
  for (int i = 0; i < D; i++) x[i] += a[i];
  return *this;
}

//==============================================================================
template <class T, int D>
Vector<T,D>& Vector<T,D>::operator-=(Vector<T,D> a)
{
  for (int i = 0; i < D; i++) x[i] -= a[i];
  return *this;
}

//==============================================================================
template <class T, int D>
Vector<T,D>& Vector<T,D>::operator*=(T l)
{
  for (int i = 0; i < D; i++) x[i] *= l;
  return *this;
}

//==============================================================================
template <class T, int D>
Vector<T,D> operator+(Vector<T,D> a, Vector<T,D> b)
{
  Vector<T,D> t = 0;
  return (t += a) += b;
}

//==============================================================================
template <class T, int D>
Vector<T,D> operator-(Vector<T,D> a)
{
  Vector<T,D> t = 0;
  return t -= a;
}

//==============================================================================
template <class T, int D>
Vector<T,D> operator-(Vector<T,D> a, Vector<T,D> b)
{
  Vector<T,D> t;
  return t = a + (-b);
}

//==============================================================================
template <class T, int D>
Vector<T,D> operator*(T l, Vector<T,D> a)
{
  Vector<T,D> t;
  return (t = a) *= l;
}

//==============================================================================
template <class T, int D>
double norm(Vector<T,D> a)
{
  double t = 0;
  for (int i = 0; i < D; i++) t += a[i]*a[i];
  return sqrt(t);
}

//==============================================================================
template <class T, int D>
double norm2(Vector<T,D> a)
{
  double t = 0;
  for (int i = 0; i < D; i++) t += a[i]*a[i];
  return t;
}

//==============================================================================
template <class T, int D>
ostream& operator<<(ostream& os, Vector<T,D> a)
{
  for (int i = 0; i < D; i++) os << a[i] << " ";
  return os;
}

//==============================================================================
template <class T, int D>
double inner (Vector<T,D> a, Vector<T,D> b)
{
  double t = 0;
  for (int i = 0; i < D; i++) t += a[i]*b[i];
  return t;
}

#endif