#ifndef VECTOR_CPP
#define VECTOR_CPP

#include <cmath>

#include "vector.h"

///\brief Copy assignment operator for Vector class. Accepts a Vector of type U
/// and casts it to type T. This is useful, e.g. for converting between Vector
/// types.
///\tparam U Type of Vector to copy.
///\tparam D Dimension of Vector to copy. Must be the same as this Vector.
///\tparam T Type of this Vector. Type U will be cast to this type.
///\param a Vector to copy.
///\return Reference to this Vector.
template <class T, int D>
template <class U>
Vector<T,D>& Vector<T,D>::operator=(const Vector<U,D> & a)
{
  for (int i = 0; i < D; i++) x[i] = (T)a(i);
  return *this;
}

///\brief Copy assignment operator for Vector class. Accepts a scalar of type
/// U and casts it to type T. This is useful for initializing a Vector with a
/// single value.
///\tparam U Type of scalar to assign to all elements of this Vector.
///\tparam D Dimension of this Vector.
///\tparam T Type of this Vector. Input double will be cast to this type.
///\param a Value to set all elements of this Vector to.
///\return Reference to this Vector.
template <class T, int D>
Vector<T,D>& Vector<T,D>::operator=(double a)
{
  for (int i = 0; i < D; i++) x[i] = (T)a;
  return *this;
}

///\brief Addition assignment operator for Vector class
///\tparam U Type of Vector to add to this Vector.
///\tparam D Dimension of Vector to add. Must be the same as this Vector.
///\tparam T Type of this Vector. Type U will be cast to this type.
///\param a Vector to add.
///\return Reference to this Vector.
template <class T, int D>
Vector<T,D>& Vector<T,D>::operator+=(const Vector<T,D> & a)
{
  for (int i = 0; i < D; i++) x[i] += a(i);
  return *this;
}

///\brief Subtraction assignment operator for Vector class
///\tparam U Type of Vector to subtract from this Vector.
///\tparam D Dimension of Vector to subtract. Must be the same as this Vector.
///\tparam T Type of this Vector. Type U will be cast to this type.
///\param a Vector to subtract.
///\return Reference to this Vector.
template <class T, int D>
Vector<T,D>& Vector<T,D>::operator-=(const Vector<T,D> & a)
{
  for (int i = 0; i < D; i++) x[i] -= a(i);
  return *this;
}

///\brief Multiplication assignment operator for Vector class. Multiplies
/// the entire vector for a scalar
///\tparam U Type of scalar to multiply this Vector by.
///\tparam D Dimension of this Vector.
///\tparam T Type of this Vector. Input double will be cast to this type.
///\param a Value to multiply this Vector by.
///\return Reference to this Vector.
template <class T, int D>
Vector<T,D>& Vector<T,D>::operator*=(T l)
{
  for (int i = 0; i < D; i++) x[i] *= l;
  return *this;
}

///\brief Sum operator of two Vectors.
///\tparam T Type of Vectors involved in the operation.
///\tparam D Dimension of Vectors involved in the operation.
///\param a First Vector in the sum.
///\param b Second Vector in the sum.
///\return Sum of the two Vectors.
template <class T, int D>
Vector<T,D> operator+(const Vector<T,D>& a, const Vector<T,D>& b)
{
  Vector<T,D> t;
  for (int i = 0; i < D; i++) t(i) = a(i) + b(i);
  return t;
}

///\brief Negative operator of a vector.
///\tparam T Type of the Vector involved in the operation.
///\tparam D Dimension of the  Vector involved in the operation.
///\param a vector which will be multiplied by -1.
///\return Negative of the input Vector.
template <class T, int D>
Vector<T,D> operator-(const Vector<T,D>& a)
{
  return (-1.0)*a;
}

///\brief Difference operator of two Vectors.
///\tparam T Type of Vectors involved in the operation.
///\tparam D Dimension of Vectors involved in the operation.
///\param a First Vector in the difference.
///\param b Second Vector in the difference.
///\return Difference of the two Vectors.
template <class T, int D>
Vector<T,D> operator-(const Vector<T,D>& a, const Vector<T,D>& b)
{
  Vector<T,D> t;
  for (int i = 0; i < D; i++) t(i) = a(i) - b(i);
  return t;
}

///\brief Multiplication operator of a Vector and a scalar.
///\tparam T Type of the Vector involved in the operation.
///\tparam D Dimension of the  Vector involved in the operation.
///\param a Vector to multiply.
///\param l Scalar to multiply by.
///\return Product of the Vector and the scalar.
template <class T, int D>
Vector<T,D> operator*(T l, const Vector<T,D>& a)
{
  Vector<T,D> t;
  for (int i = 0; i < D; i++) t(i) = l*a(i);
  return t;
}

///\brief Gets the norm of a Vector.
///\tparam T Type of the Vector involved in the operation.
///\tparam D Dimension of the  Vector involved in the operation.
///\param a Vector to get the norm of.
///\return Norm of the Vector.
template <class T, int D>
double Norm(const Vector<T,D>& a)
{
  double t = 0;
  for (int i = 0; i < D; i++) t += a(i)*a(i);
  return sqrt(t);
}

/// \brief Gets the squared norm of a Vector.
/// \tparam T Type of the Vector involved in the operation.
/// \tparam D Dimension of the  Vector involved in the operation.
/// \param a Vector to get the norm of.
/// \return Norm of the Vector.
template <class T, int D>
double Norm2(const Vector<T,D>& a)
{
  double t = 0;
  for (int i = 0; i < D; i++) t += a(i)*a(i);
  return t;
}

///\brief Insertion operator for Vector class.
///\tparam T Type of the Vector involved in the operation.
///\tparam D Dimension of the  Vector involved in the operation.
///\param os Output stream to insert into.
///\param a Vector to insert into the output stream.
template <class T, int D>
ostream& operator<<(ostream& os, const Vector<T,D>& a)
{
  for (int i = 0; i < D; i++) os << a(i) << " ";
  return os;
}

///\brief Inner product of two Vectors.
///\tparam T Type of Vectors involved in the operation.
///\tparam D Dimension of Vectors involved in the operation.
///\param a First Vector in the inner product.
///\param b Second Vector in the inner product.
///\return Inner product of the two Vectors.
template <class T, int D>
double inner (const Vector<T,D>& a, const Vector<T,D>& b)
{
  double t = 0;
  for (int i = 0; i < D; i++) t += a(i)*b(i);
  return t;
}

#endif