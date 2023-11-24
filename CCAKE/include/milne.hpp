#ifndef MILNE_H
#define MILNE_H

#include <Kokkos_Core.hpp>

namespace ccake{
namespace milne{

template <class T, int D>
class Vector
{
private:
  T x[D];
  int status = 1; ///< Stores if the vector is covariant or contravariant. Contravariant by default.
                  ///TODO: I think this may not be necessary
public:
  KOKKOS_FUNCTION Vector<T,D>(){for (int i = 0; i < D; i++) x[i] = (T)0.0;}
  KOKKOS_FUNCTION Vector<T,D>(T x0){for (int i = 0; i < D; i++) x[i] = x0;}
  KOKKOS_FUNCTION Vector<T,D>(std::initializer_list<T> x0){
    #ifdef DEBUG
    assert(x0.size() == D);
    #endif
    std::copy(x0.begin(), x0.end(), x);
  }

  KOKKOS_FUNCTION T& operator()(const int i)       { return x[i]; };
  KOKKOS_FUNCTION T  operator()(const int i) const { return x[i]; };
  KOKKOS_FUNCTION Vector<T,D>& operator=(const Vector<T,D>&);
  template <class U> KOKKOS_FUNCTION
  Vector<T,D>& operator=(U);
  template <class U> KOKKOS_FUNCTION
  Vector<T,D>& operator+=(const Vector<U,D>&);
  template <class U> KOKKOS_FUNCTION
  Vector<T,D>& operator-=(const Vector<U,D>&);
  template <class U> KOKKOS_FUNCTION
  Vector<T,D>& operator*=(U);
  KOKKOS_FUNCTION void make_covariant(double t2) {x[D-1] *= t2;status-=2;};
  KOKKOS_FUNCTION void make_contravariant(double t2){x[D-1] /= t2;status+=2;};
  KOKKOS_FUNCTION void set_covariance(int c){ status = c;}; //Adjust the covariance of the object
  KOKKOS_FUNCTION ~Vector<T,D>(){};
};

//==============================================================================
// overloaded operations +, - and * between scalars and vectors
template <class T, int D> KOKKOS_FUNCTION 
Vector<T,D> operator+ ( const Vector<T,D>&, const Vector<T,D>& );
template <class T, int D> KOKKOS_FUNCTION 
Vector<T,D> operator- ( const Vector<T,D>&, const Vector<T,D>& );
template <class T, int D> KOKKOS_FUNCTION 
Vector<T,D> operator* ( T, const Vector<T,D>&           );


template <class T, int D1, int D2>
class Matrix{
/// @class Stores data as a D1 by D2 matrix
/// @details We assume that D1 are the number of rows and D2 the number
/// of columns. When we invoking D(i,j), we are retrieving the element
/// in the i-th row and j-th column. The fast running index is the last 
/// one (the j index)
private:
  T x[D1*D2];
  int status[2] = {1,1}; ///< Stores if the matrix indices are covariant or contravariant. All contravariant by default.
                         ///TODO: I think this may not be necessary
  KOKKOS_INLINE_FUNCTION int index(const int i, const int j) const { return i*D2+j; };
public:
  KOKKOS_FUNCTION Matrix<T,D1,D2>() {for (int i=0; i<D1*D2; i++) x[i] = 0;};
  KOKKOS_FUNCTION Matrix<T,D1,D2>(T x0) {for (int i = 0; i < D1*D2; i++) x[i] = x0;};
  KOKKOS_FUNCTION T operator()(const int i, const int j) const { return x[index(i,j)]; };
  KOKKOS_FUNCTION T& operator()(const int i, const int j) { return x[index(i,j)];};
  KOKKOS_FUNCTION Matrix<T,D1,D2>& operator=(T); //Documented
  template <class U>
  KOKKOS_FUNCTION Matrix<T,D1,D2>& operator=(const Matrix<U,D1,D2>& a); //< Implemented
  KOKKOS_FUNCTION Matrix<T,D1,D2>& operator+=(const Matrix<T,D1,D2>&);
  KOKKOS_FUNCTION Matrix<T,D1,D2>& operator-=(const Matrix<T,D1,D2>&);
  KOKKOS_FUNCTION Matrix<T,D1,D2>& operator*=(T);
  KOKKOS_FUNCTION void make_covariant(int slot, double t2);
  KOKKOS_FUNCTION void make_contravariant(int slot, double t2);
  KOKKOS_FUNCTION void set_covariance(int c[2]){ status[0] = c[0]; status[1] = c[1];}; //Adjust the covariance of the object
  KOKKOS_FUNCTION ~Matrix<T,D1,D2>(){};
};

//----------------------------------
//Overload operators +, - and * between scalars, vectors and matrices

//==============================
// sum two matrices
template <class T, int D1, int D2> KOKKOS_FUNCTION
Matrix<T,D1,D2> operator+ ( const Matrix<T,D1,D2>&, const Matrix<T,D1,D2>& );
//==============================
// negative of Matrix
template <class T, int D1, int D2> KOKKOS_FUNCTION
Matrix<T,D1,D2> operator- ( const Matrix<T,D1,D2>& );
//==============================
// subtract two matrices
template <class T, int D1, int D2> KOKKOS_FUNCTION
Matrix<T,D1,D2> operator- ( const Matrix<T,D1,D2>&, const Matrix<T,D1,D2>& );
//==============================
// multiply Matrix by scalar
template <class T, int D1, int D2> KOKKOS_FUNCTION
Matrix<T,D1,D2> operator* ( T, const Matrix<T,D1,D2>& );
//==============================
// standard Matrix multiplication
template <class T, int D1, int D2, int Da2, int Db1> KOKKOS_FUNCTION
Matrix<T,D1,D2> operator* ( const Matrix<T,D1,Da2>& a, const Matrix<T,Db1,D2>& b );
//==============================
// right-multiply Matrix a on vector b: a.b
template <class T, int D1, int D2> KOKKOS_FUNCTION
Vector<T,D1>    operator* ( const Matrix<T,D1,D2>& a, const Vector<T,D2>& b );
//==============================
// left-multiply Matrix a on vector b: b.a
template <class T, int D1, int D2> KOKKOS_FUNCTION
Vector<T,D2>    operator* ( const Vector<T,D1>& a, const Matrix<T,D1,D2>& b );

//==============================================================================
// Auxiliary operations
template <class T, int D> KOKKOS_FUNCTION
double inner( const Vector<T,D>&, const Vector<T,D>& );
template <class T, int D1, int D2> KOKKOS_FUNCTION
Matrix<T,D2,D1> transpose( const Matrix<T, D1, D2>& a );
template <class T, int D> KOKKOS_FUNCTION
Matrix<T,D,D> inverse(const Matrix<T, D, D>& a);
template <class T, int D1, int D2> KOKKOS_FUNCTION
Vector<T,(D1-1)> colp1( int l, const Matrix<T, D1, D2>& a );
template <class T, int D1, int D2> KOKKOS_FUNCTION
void mini(Matrix<T,D1-1,D2-1> &b, const Matrix<T,D1,D2>& a );


//================================================================
//Implementations of the Vector class 
//================================================================

///\brief Copy assignment operator for Vector class. Accepts a Vector of type U
/// and casts it to type T. This is useful, e.g. for converting between Vector
/// types.
///\tparam U Type of Vector to copy.
///\tparam D Dimension of Vector to copy. Must be the same as this Vector.
///\tparam T Type of this Vector. Type U will be cast to this type.
///\param a Vector to copy.
///\return Reference to this Vector.
template <class T, int D> KOKKOS_FUNCTION
Vector<T,D>& Vector<T,D>::operator=(const Vector<T,D> & a)
{
  for (int i = 0; i < D; i++) x[i] = a(i);
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
template <class T, int D> template <class U> KOKKOS_FUNCTION
Vector<T,D>& Vector<T,D>::operator=(U a)
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
template <class T, int D> template <class U> KOKKOS_FUNCTION
Vector<T,D>& Vector<T,D>::operator+=(const Vector<U,D> & a)
{
  for (int i = 0; i < D; i++) x[i] += (T)a(i);
  return *this;
}

///\brief Subtraction assignment operator for Vector class
///\tparam U Type of Vector to subtract from this Vector.
///\tparam D Dimension of Vector to subtract. Must be the same as this Vector.
///\tparam T Type of this Vector. Type U will be cast to this type.
///\param a Vector to subtract.
///\return Reference to this Vector.
template <class T, int D> template <class U> KOKKOS_FUNCTION
Vector<T,D>& Vector<T,D>::operator-=(const Vector<U,D> & a)
{
  for (int i = 0; i < D; i++) x[i] -= (T)a(i);
  return *this;
}

///\brief Multiplication assignment operator for Vector class. Multiplies
/// the entire vector for a scalar
///\tparam U Type of scalar to multiply this Vector by.
///\tparam D Dimension of this Vector.
///\tparam T Type of this Vector. Input double will be cast to this type.
///\param a Value to multiply this Vector by.
///\return Reference to this Vector.
template <class T, int D> template <class U> KOKKOS_FUNCTION
Vector<T,D>& Vector<T,D>::operator*=(U a)
{
  for (int i = 0; i < D; i++) x[i] *= (T)a;
  return *this;
}

//Dimension 2 is special in the sense we do not distinguish
//between covariant and contravariant components
template <> KOKKOS_INLINE_FUNCTION
void Vector<double,2>::make_covariant(double t2) {status-=2;};
template <> KOKKOS_INLINE_FUNCTION
void Vector<float,2>::make_covariant(double t2) {status-=2;};
template <> KOKKOS_INLINE_FUNCTION
void Vector<int,2>::make_covariant(double t2) {status-=2;};
template <> KOKKOS_INLINE_FUNCTION
void Vector<double,2>::make_contravariant(double t2) {status+=2;};
template <> KOKKOS_INLINE_FUNCTION
void Vector<float,2>::make_contravariant(double t2) {status+=2;};
template <> KOKKOS_INLINE_FUNCTION
void Vector<int,2>::make_contravariant(double t2) {status+=2;};

//================================================================
//Implementations of the Matrix class 
//================================================================

/// @brief Assign a scalar to every entry of the matrix
/// @tparam T the data type of the matrix
/// @tparam D1 The number of rows in the matrix
/// @tparam D2 The number of columns in the matrix
/// @param a A d
/// @return The pointer to the copied matrix 
template <class T, int D1, int D2> KOKKOS_FUNCTION
Matrix<T,D1,D2>& Matrix<T,D1,D2>::operator=(T a)
{
  int k = 0;
  for (int i=0; i<D1; i++)
  for (int j=0; j<D2; j++)
    x[k++] = (T)a;
  return *this; // AOK
}

/// @brief Copy matrix A into matrix B
/// @tparam T the data type of the matrix
/// @tparam D1 The number of rows in the matrix
/// @tparam D2 The number of columns in the matrix
template <class T, int D1, int D2> template<class U> KOKKOS_FUNCTION
Matrix<T,D1,D2>& Matrix<T,D1,D2>::operator=(const Matrix<U,D1,D2>& a)
{
  int k = 0;
  for (int i=0; i<D1; i++)
  for (int j=0; j<D2; j++)
    x[k++] = (T)a(i,j);
  return *this; // AOK
}


template <class T, int D1, int D2> KOKKOS_FUNCTION
Matrix<T,D1,D2>& Matrix<T,D1,D2>::operator+=(const Matrix<T,D1,D2>& a)
{
  int k = 0;
  for(int i=0; i<D1; i++)
  for(int j=0; j<D2; j++)
    x[k++] += a(i,j);
  return *this; // AOK
}


template <class T, int D1, int D2> KOKKOS_FUNCTION
Matrix<T,D1,D2>& Matrix<T,D1,D2>::operator-=(const Matrix<T,D1,D2>& a)
{
  int k = 0;
  for(int i=0; i<D1; i++)
  for(int j=0; j<D2; j++)
    x[k++] -= a(i,j);
  return *this; // AOK
}


template <class T, int D1, int D2> KOKKOS_FUNCTION
Matrix<T,D1,D2>& Matrix<T,D1,D2>::operator*=(T l)
{
  for (int i=0; i<D1*D2; i++) x[i] *= l;
  return *this; // AOK
}

/// @brief Lower one index of the matrix
/// @tparam T the data type of the matrix
/// @tparam D1 The number of rows in the matrix
/// @tparam D2 The number of columns in the matrix
/// @param slot Which index to lower.
/// @param The square of the temporal coordinate: (x^0)^2
/// @details slot values can be 0 or 1. Zero means we are lowering
/// the first index, i.e. the row index. 1 is lowering the second index
/// (the column index). Notice the absence of minus sign. This is because
/// we assume the metric does not couple space and time components and
/// we need to operate only in the space components
template <class T, int D1, int D2> KOKKOS_FUNCTION
void Matrix<T,D1,D2>::make_covariant(int slot, double t2){
  status[slot]-=2;
  switch (slot){
    case 0: ///Lower first index. Array becomes M_i^{\,j} = g_ik M^{kj} = tau^2 M^{3j}\delta_i^3
      for (int i=0; i<D1; ++i) x[index(D1-1,i)] *= t2;
      break;
    case 1:
      for (int i=0; i<D2; ++i) x[index(i,D2-1)] *= t2;
      break;
    default:
      break;
  }
}
/// @brief Raise one index of the matrix
/// @tparam T the data type of the matrix
/// @tparam D1 The number of rows in the matrix
/// @tparam D2 The number of columns in the matrix
/// @param slot Which index to lower.
/// @param The square of the temporal coordinate: (x^0)^2
/// @details slot values can be 0 or 1. Zero means we are lowering
/// the first index, i.e. the row index. 1 is lowering the second index
/// (the column index). Notice the absence of minus sign. This is because
/// we assume the metric does not couple space and time components and
/// we need to operate only in the space components
template <class T, int D1, int D2> KOKKOS_FUNCTION
void Matrix<T,D1,D2>::make_contravariant(int slot, double t2){
  status[slot]+=2;
  switch (slot){
    case 0: 
      for (int i=0; i<D1; ++i) x[index(D1-1,i)] /= t2;
      break;
    case 1:
      for (int i=0; i<D2; ++i) x[index(i,D2-1)] /= t2;
      break;
    default:
      break;
  }
}

/// 2x2 matrices are special. We do not distinguish covariant and contravariant ones

template <> KOKKOS_INLINE_FUNCTION
void Matrix<double,2,2>::make_contravariant(int slot, double t2){
  status[slot]+=2;
}
template <> KOKKOS_INLINE_FUNCTION
void Matrix<float,2,2>::make_contravariant(int slot, double t2){
  status[slot]+=2;
}
template <> KOKKOS_INLINE_FUNCTION
void Matrix<int,2,2>::make_contravariant(int slot, double t2){
  status[slot]+=2;
}
template <> KOKKOS_INLINE_FUNCTION
void Matrix<double,2,2>::make_covariant(int slot, double t2){
  status[slot]-=2;
}
template <> KOKKOS_INLINE_FUNCTION
void Matrix<float,2,2>::make_covariant(int slot, double t2){
  status[slot]-=2;
}
template <> KOKKOS_INLINE_FUNCTION
void Matrix<int,2,2>::make_covariant(int slot, double t2){
  status[slot]-=2;
}

//================================================================
//Implementations of vector overload operations
//================================================================
///\brief Sum operator of two Vectors.
///\tparam T Type of Vectors involved in the operation.
///\tparam D Dimension of Vectors involved in the operation.
///\param a First Vector in the sum.
///\param b Second Vector in the sum.
///\return Sum of the two Vectors.
template <class T, int D> KOKKOS_FUNCTION
Vector<T,D> operator+(const Vector<T,D>& a, const Vector<T,D>& b)
{
  Vector<T,D> t;
  for (int i = 0; i < D; i++) t(i) = a(i) + b(i);
  return t;
}

///\brief Difference operator of two Vectors.
///\tparam T Type of Vectors involved in the operation.
///\tparam D Dimension of Vectors involved in the operation.
///\param a First Vector in the difference.
///\param b Second Vector in the difference.
///\return Difference of the two Vectors.
template <class T, int D> KOKKOS_FUNCTION
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
template <class T, int D> KOKKOS_FUNCTION
Vector<T,D> operator*(T l, const Vector<T,D>& a)
{
  Vector<T,D> t;
  for (int i = 0; i < D; i++) t(i) = l*a(i);
  return t;
}

//================================================================
//Implementations of Matrix overload operations
//================================================================

template <class T, int D1, int D2> KOKKOS_FUNCTION
Matrix<T,D1,D2> operator+(const Matrix<T,D1,D2>& a, const Matrix<T,D1,D2>& b)
{
  Matrix<T,D1,D2> t;
  for (int i=0; i<D1; i++)
  for (int j=0; j<D2; j++)
    t(i,j) = a(i,j)+b(i,j);
  return t;
} // AOK


template <class T, int D1, int D2> KOKKOS_FUNCTION
Matrix<T,D1,D2> operator-(const Matrix<T,D1,D2>& a)
{
  return (-1.0)*a;
} // AOK


template <class T, int D1, int D2> KOKKOS_FUNCTION
Matrix<T,D1,D2> operator-(const Matrix<T,D1,D2>& a, const Matrix<T,D1,D2>& b)
{
  Matrix<T,D1,D2> t;
  for (int i = 0; i < D1; i++)
  for (int j = 0; j < D2; j++)
    t(i,j) = a(i,j) - b(i,j);
  return t;
} // AOK

/// @brief  Standard matrix multiplication
/// @tparam T Data type of the matrices
/// @tparam D1 Number of rows in matrix a
/// @tparam D2 Number of columns in matrix b
/// @tparam Dk Number of columns in matrix a and rows in matrix b
/// @param a Input matrix a
/// @param b Input matrix b
/// @details Perform the multiplication C(i,k) = Sum_{k=1,Dk}a(i,k)b(k,j)
/// @return The matrix C(i,k)
template <class T, int D1, int D2, int Dk> KOKKOS_FUNCTION
Matrix<T,D1,D2> operator*(const Matrix<T,D1,Dk>& a, const Matrix<T,Dk,D2>& b)
{

  Matrix<T,D1,D2> t;
  for (int i=0; i<D1; i++)
  for (int j=0; j<D2; j++)
  for (int k=0; k<Dk; k++)
    t(i,j) += a(i,k) * b(k,j);

  return t;
} // AOK


template <class T, int D1, int D2> KOKKOS_FUNCTION
Matrix<T,D1,D2> operator*(T l, const Matrix<T,D1,D2>& a)
{
  Matrix<T,D1,D2> t;
  for (int i = 0; i < D1; i++)
  for (int j = 0; j < D2; j++)
    t(i,j) = l*a(i,j);

  return t;
} // AOK

/// @brief 
/// @tparam T 
/// @tparam D1 
/// @tparam D2 
/// @param a 
/// @param b 
/// @return 
template <class T, int D1, int D2> KOKKOS_FUNCTION
Vector<T,D1> operator*(const Matrix<T,D1,D2>& a, const Vector<T,D2>& b)
{
  Vector<T,D1> t;
  for (int i=0; i<D1; i++)
  for (int j=0; j<D2; j++)
    t(i) += a(i,j) * b(j);

  return t;
} // AOK

template <class T, int D1, int D2> KOKKOS_FUNCTION
Vector<T,D2> operator*(const Vector<T,D1>& b, const Matrix<T,D1,D2>& a)
{
  Vector<T,D1> t;
  for (int j=0; j<D2; j++)
  for (int i=0; i<D1; i++)
    t(j) += b(i)*a(i,j) ;

  return t;
} // AOK



template <class T, int D1, int D2> KOKKOS_FUNCTION
Matrix<T,D1,D2> operator* ( const Vector<T,D1>& a, const Vector<T,D2>& b )
{
  Matrix<T,D1,D2> t;
  for (int i = 0; i < D1; i++)
  for (int j = 0; j < D2; j++)
    t(i,j) = a(i)*b(j);
  return t;
  
} // AOK


//================================================================
//Implementations of auxiliary operations
//================================================================
///\brief Inner product of two Vectors.
///\tparam T Type of Vectors involved in the operation.
///\tparam D Dimension of Vectors involved in the operation.
///\param a First Vector in the inner product.
///\param b Second Vector in the inner product.
///\return Inner product of the two Vectors.
template <class T, int D> KOKKOS_FUNCTION
double inner (const Vector<T,D>& a, const Vector<T,D>& b)
{
  double t = 0;
  for (int i = 0; i < D; i++) t += a(i)*b(i);
  return t;
}

template <class T, int D1, int D2> KOKKOS_FUNCTION
Matrix<T,D2,D1> transpose(const Matrix<T,D1,D2>& a)
{
  Matrix<T,D2,D1> t;

  for (int i=0; i<D2; i++)
  for (int j=0; j<D1; j++)
    t(i,j) = a(j,i);

  return t;
} // AOK

template<> KOKKOS_INLINE_FUNCTION
Matrix<double,1,1> inverse<double, 1>(const Matrix<double,1,1> &A)
{
  Matrix<double,1,1> Ainv;
  Ainv(0,0) = 1./A(0,0);
  return Ainv;
}

template<>
KOKKOS_INLINE_FUNCTION
Matrix<double,2,2> inverse<double, 2>(const Matrix<double,2,2> &A)
{
  Matrix<double,2,2> Ainv;
  double det = A(0,0)*A(1,1) - A(0,1)*A(1,0);
  Ainv(0,0) = A(1,1)/det;
  Ainv(0,1) = -A(0,1)/det;
  Ainv(1,0) = -A(1,0)/det;
  Ainv(1,1) = A(0,0)/det;
  return Ainv;
}

template<>
KOKKOS_INLINE_FUNCTION
Matrix<double,3,3> inverse<double, 3>(const Matrix<double,3,3> &A)
{
    Matrix<double,3,3> Ainv;
    double det = A(0,0)*A(1,1)*A(2,2) + A(0,1)*A(1,2)*A(2,0) +
                 A(0,2)*A(1,0)*A(2,1) - A(0,2)*A(1,1)*A(2,0) -
                 A(0,1)*A(1,0)*A(2,2) - A(0,0)*A(1,2)*A(2,1);
    Ainv(0,0) = (A(1,1)*A(2,2) - A(1,2)*A(2,1))/det;
    Ainv(0,1) = (A(0,2)*A(2,1) - A(0,1)*A(2,2))/det;
    Ainv(0,2) = (A(0,1)*A(1,2) - A(0,2)*A(1,1))/det;
    Ainv(1,0) = (A(1,2)*A(2,0) - A(1,0)*A(2,2))/det;
    Ainv(1,1) = (A(0,0)*A(2,2) - A(0,2)*A(2,0))/det;
    Ainv(1,2) = (A(0,2)*A(1,0) - A(0,0)*A(1,2))/det;
    Ainv(2,0) = (A(1,0)*A(2,1) - A(1,1)*A(2,0))/det;
    Ainv(2,1) = (A(0,1)*A(2,0) - A(0,0)*A(2,1))/det;
    Ainv(2,2) = (A(0,0)*A(1,1) - A(0,1)*A(1,0))/det;
    return Ainv;
}

// takes spatial components of l^{th} column in space-time tensor (Matrix a)
template <class T, int D1, int D2> KOKKOS_FUNCTION
Vector<T,(D1-1)> colp1(int l, const Matrix<T,D1,D2>& a)
{
  Vector<T,(D1-1)> v;
  for(int i=1; i<D1; i++) v(i-1)=(T)a(i,l);
  return v;
} // AOK

// store transverse spatial components of a in b
template <class T, int D1, int D2> KOKKOS_FUNCTION
void mini( Matrix<T,D1-1,D2-1> &b, const Matrix<T,D1,D2>& a)
{
  for(int i=1; i<D1; i++)
  for(int j=1; j<D2; j++)
    b(i-1,j-1) = (T)a(i,j);
} // AOK

template <class T, int D1, int D2> KOKKOS_FUNCTION
double con( const Matrix<T,D1,D2>& a, const Matrix<T,D1,D2>& b )
{
  double t = 0.0;
  for (int i = 0; i < D2; i++)
  for (int j = 0; j < D1; j++)
    t += a(j,i)*b(j,i);
  return t;
} // AOK

// takes spatial components of l^{th} row in space-time tensor (Matrix a)
template <class T, int D1, int D2> KOKKOS_FUNCTION
Vector<T,(D2-1)> rowp1(int l, const Matrix<T,D1,D2>& a)
{
  Vector<T,(D2-1)> v;
  for (int i=1; i<D2; i++) v(i-1)=(T)a(l,i);
  return v;
} // AOK

// store a in transverse spatial components of b
template <class T, int D1, int D2> KOKKOS_FUNCTION
void tmini( Matrix<T,D1,D2>& b, const Matrix<T,D1-1,D2-1>& a )
{
  for (int j = 0; j < (D2-1); j++)
  for (int i = 0; i < (D1-1); i++)
    b(i+1,j+1) = (T)a(i,j);
} // AOK

template <class T, int D> KOKKOS_INLINE_FUNCTION
T tr( const Matrix<T,D,D>& a,T t2)
{
  double s = 0.0;
  for (int i = 0; i < D-1; i++)
    s += a(i,i);
  s += a(D-1,D-1)*t2;
  exit(8);
  return s;
} // AOK

template <> KOKKOS_INLINE_FUNCTION
double tr( const Matrix<double,2,2>& a,double t2)
{
  double s = 0.0;
  for (int i = 0; i < 2; i++)
    s += a(i,i);
  return s;
} // AOK


template <class T, int D1, int D2> KOKKOS_FUNCTION
double con2(const Matrix<T,D1,D2>& a, const Matrix<T,D1,D2>& b)
{
  double t = 0.0;
  for (int i = 0; i < D2; i++)
  for (int j = 0; j < D1; j++)
    t += a(i,j)*b(i,j);
  return t;
} // AOK

}}
#endif