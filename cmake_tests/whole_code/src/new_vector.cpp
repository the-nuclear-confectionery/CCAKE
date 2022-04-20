#ifndef NEW_VECTOR_CPP
#define NEW_VECTOR_CPP

#include <cmath>

#include "new_vector.h"

//==============================================================================
template <class T, int D>
template <class U>
Vector<T,D>& Vector<T,D>::operator=(const Vector<U,D> & a)
{
  for (int i = 0; i < D; i++) x[i] = (T)a(i);
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
Vector<T,D>& Vector<T,D>::operator+=(const Vector<T,D> & a)
{
  for (int i = 0; i < D; i++) x[i] += a(i);
  return *this;
}

//==============================================================================
template <class T, int D>
Vector<T,D>& Vector<T,D>::operator-=(const Vector<T,D> & a)
{
  for (int i = 0; i < D; i++) x[i] -= a(i);
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
Vector<T,D> operator+(const Vector<T,D>& a, const Vector<T,D>& b)
{
  Vector<T,D> t;
  for (int i = 0; i < D; i++) t(i) = a(i) + b(i);
  return t;
}

//==============================================================================
template <class T, int D>
Vector<T,D> operator-(const Vector<T,D>& a)
{
  return (-1.0)*a;
}

//==============================================================================
template <class T, int D>
Vector<T,D> operator-(const Vector<T,D>& a, const Vector<T,D>& b)
{
  Vector<T,D> t;
  for (int i = 0; i < D; i++) t(i) = a(i) - b(i);
  return t;
}

//==============================================================================
template <class T, int D>
Vector<T,D> operator*(T l, const Vector<T,D>& a)
{
  Vector<T,D> t;
  for (int i = 0; i < D; i++) t(i) = l*a(i);
  return t;
}

//==============================================================================
template <class T, int D>
double norm(const Vector<T,D>& a)
{
  double t = 0;
  for (int i = 0; i < D; i++) t += a(i)*a(i);
  return sqrt(t);
}

//==============================================================================
template <class T, int D>
double norm2(const Vector<T,D>& a)
{
  double t = 0;
  for (int i = 0; i < D; i++) t += a(i)*a(i);
  return t;
}

//==============================================================================
template <class T, int D>
ostream& operator<<(ostream& os, const Vector<T,D>& a)
{
  for (int i = 0; i < D; i++) os << a(i) << " ";
  return os;
}

//==============================================================================
template <class T, int D>
double inner (const Vector<T,D>& a, const Vector<T,D>& b)
{
  double t = 0;
  for (int i = 0; i < D; i++) t += a(i)*b(i);
  return t;
}

#endif