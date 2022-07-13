#include <iomanip>
#include <iostream>

#include "colors.h"
#include "matrix.h"
#include "vector.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline void passed(const string & f) { cout << setw(50) << left << " "+f << ": " << GREEN << "PASSED\n" << RESET; }
inline void failed(const string & f) { cout << setw(50) << left << " "+f << ": " << RED << "FAILED\n" << RESET; }

template <class T, int D>
bool operator==( const Vector<T,D>& a, const Vector<T,D>& b )
{
  bool result = true;
  for (int i = 0; i < D; i++) result = result && ( a(i) == b(i) );
  return result;
}


template <class T, int D1, int D2>
bool operator==( const Matrix<T,D1,D2>& a, const Matrix<T,D1,D2>& b )
{
  bool result = true;
  for (int i = 0; i < D1; i++)
  for (int j = 0; j < D2; j++)
    result = result && ( a(i,j) == b(i,j) );
  return result;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
namespace vector_checks
{
  ////////////////////////////////////////////////////////////////////////////////
  void check_copy()
  {
    Vector<double, 2> a, b;
    a(0) = 1; a(1) = 2;
    b = a;
    if ( a == b ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }

  ////////////////////////////////////////////////////////////////////////////////
  void check_assignment()
  {
    Vector<double, 2> a = 7;
    if ( a(0) == 7 && a(1) == 7 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }



  ////////////////////////////////////////////////////////////////////////////////
  void check_inplace_addition()
  {
    Vector<double, 2> a = 7., b = 1.;
    a += b;
    if ( a(0) == 8 && a(1) == 8 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }

  ////////////////////////////////////////////////////////////////////////////////
  void check_inplace_subtraction()
  {
    Vector<double, 2> a = 7., b = 1.;
    a -= b;
    if ( a(0) == 6 && a(1) == 6 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }


  ////////////////////////////////////////////////////////////////////////////////
  void check_inplace_multiplication()
  {
    Vector<double, 2> a = 7.;
    double b = 3;
    a *= b;
    if ( a(0) == 21 && a(1) == 21 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }

  ////////////////////////////////////////////////////////////////////////////////
  void check_negation()
  {
    Vector<double, 2> a = 7.;
    a = -a;
    if ( a(0) == -7 && a(1) == -7 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }


  ////////////////////////////////////////////////////////////////////////////////
  void check_addition()
  {
    Vector<double, 2> a = 7, b = 3;
    Vector<double, 2> c = a + b;
    if ( c(0) == 10 && c(1) == 10 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }


  ////////////////////////////////////////////////////////////////////////////////
  void check_subtraction()
  {
    Vector<double, 2> a = 7, b = 3;
    Vector<double, 2> c = a - b;
    if ( c(0) == 4 && c(1) == 4 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }


  ////////////////////////////////////////////////////////////////////////////////
  void check_multiplication()
  {
    Vector<double, 2> a = 7;
    double b = 3;
    Vector<double, 2> c = b*a;
    if ( c(0) == 21 && c(1) == 21 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }



  ////////////////////////////////////////////////////////////////////////////////
  void check_inner()
  {
    Vector<double, 2> a = 7;
    Vector<double, 2> b = 3;
    a(1) -= 1; b(0) += 1;
    if ( inner(a,b) == 46 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }


  ////////////////////////////////////////////////////////////////////////////////
  void check_Norm()
  {
    Vector<double, 2> a = 3;
    a(1) += 1;
    if ( Norm(a) == 5 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }




  ////////////////////////////////////////////////////////////////////////////////
  void check_Norm2()
  {
    Vector<double, 2> a = 3;
    a(1) += 1;
    if ( Norm2(a) == 25 ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }

}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void run_vector_tests()
{
  cout << "Running vector checks" << endl;
  using namespace vector_checks;

  check_copy();
  check_assignment();
  check_inplace_addition();
  check_inplace_subtraction();
  check_inplace_multiplication();
  check_negation();
  check_addition();
  check_subtraction();
  check_multiplication();
  check_inner();
  check_Norm();
  check_Norm2();

  cout << endl;
}






////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
namespace matrix_checks
{
  ////////////////////////////////////////////////////////////////////////////////
  void check_copy()
  {
    Matrix<double, 2, 2> a, b;
    a(0,0) = 1; a(0,1) = 2; a(1,0) = 3; a(1,1) = 4;
    b = a;
    if ( a == b ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }

  ////////////////////////////////////////////////////////////////////////////////
  void check_assignment()
  {
    Matrix<double, 2, 2> a = 7;
    bool test_condition = (a(0,0) == 7 && a(0,1) == 7 && a(1,0) == 7 && a(1,1) == 7);
    if ( test_condition ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }



  ////////////////////////////////////////////////////////////////////////////////
  void check_inplace_addition()
  {
    Matrix<double, 2, 2> a = 7., b = 1.;
    a += b;
    bool test_condition = (a(0,0) == 8 && a(0,1) == 8 && a(1,0) == 8 && a(1,1) == 8);
    if ( test_condition ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }

  ////////////////////////////////////////////////////////////////////////////////
  void check_inplace_subtraction()
  {
    Matrix<double, 2, 2> a = 7., b = 1.;
    a -= b;
    bool test_condition = (a(0,0) == 6 && a(0,1) == 6 && a(1,0) == 6 && a(1,1) == 6);
    if ( test_condition ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }


  ////////////////////////////////////////////////////////////////////////////////
  void check_inplace_scalar_multiplication()
  {
    Matrix<double, 2, 2> a = 7.;
    double b = 3;
    a *= b;
    bool test_condition = (a(0,0) == 21 && a(0,1) == 21 && a(1,0) == 21 && a(1,1) == 21);
    if ( test_condition ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }

  ////////////////////////////////////////////////////////////////////////////////
  void check_negation()
  {
    Matrix<double, 2, 2> a = 7.;
    a = -a;
    bool test_condition = (a(0,0) == -7 && a(0,1) == -7 && a(1,0) == -7 && a(1,1) == -7);
    if ( test_condition ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }


  ////////////////////////////////////////////////////////////////////////////////
  void check_addition()
  {
    Matrix<double, 2, 2> a = 7, b = 3;
    Matrix<double, 2, 2> c = a + b;
    bool test_condition = (c(0,0) == 10 && c(0,1) == 10 && c(1,0) == 10 && c(1,1) == 10);
    if ( test_condition ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }


  ////////////////////////////////////////////////////////////////////////////////
  void check_subtraction()
  {
    Matrix<double, 2, 2> a = 7, b = 3;
    Matrix<double, 2, 2> c = a - b;
    bool test_condition = (c(0,0) == 4 && c(0,1) == 4 && c(1,0) == 4 && c(1,1) == 4);
    if ( test_condition ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }


  ////////////////////////////////////////////////////////////////////////////////
  void check_scalar_multiplication()
  {
    Matrix<double, 2, 2> a = 7;
    double b = 3;
    Matrix<double, 2, 2> c = b*a;
    bool test_condition = (c(0,0) == 21 && c(0,1) == 21 && c(1,0) == 21 && c(1,1) == 21);
    if ( test_condition ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }

  ////////////////////////////////////////////////////////////////////////////////
  void check_matrix_matrix_multiplication()
  {
    Matrix<double, 2, 2> a;
    a(0,0) = 1; a(0,1) = 2; a(1,0) = 3; a(1,1) = 4;
    Matrix<double, 2, 2> b = 5;
    b += a;
    Matrix<double, 2, 2> c = a*b;
    bool test_condition = (c(0,0) == 22 && c(0,1) == 25 && c(1,0) == 50 && c(1,1) == 57);
    if ( test_condition ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }

  ////////////////////////////////////////////////////////////////////////////////
  void check_vector_matrix_multiplication()
  {
    Matrix<double, 2, 2> a;
    a(0,0) = 1; a(0,1) = 2; a(1,0) = 3; a(1,1) = 4;
    Vector<double, 2> b = 5; b(1) = 6;
    Vector<double, 2> c = b*a;
    bool test_condition = (c(0) == 23 && c(1) == 34);
    if ( test_condition ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }

  ////////////////////////////////////////////////////////////////////////////////
  void check_matrix_vector_multiplication()
  {
    Matrix<double, 2, 2> a;
    a(0,0) = 1; a(0,1) = 2; a(1,0) = 3; a(1,1) = 4;
    Vector<double, 2> b = 5; b(1) = 6;
    Vector<double, 2> c = a*b;
    bool test_condition = (c(0) == 17 && c(1) == 39);
    if ( test_condition ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }

  ////////////////////////////////////////////////////////////////////////////////
  void check_transpose()
  {
    Matrix<double, 2, 2> a;
    a(0,0) = 1; a(0,1) = 2; a(1,0) = 3; a(1,1) = 4;
    Matrix<double, 2, 2> aT;
    aT(0,0) = 1; aT(0,1) = 3; aT(1,0) = 2; aT(1,1) = 4;

    bool test_condition = (aT == transpose(a));
    if ( test_condition ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }

  ////////////////////////////////////////////////////////////////////////////////
  void check_determinant_2D()
  {
    Matrix<double, 2, 2> a;
    a(0,0) = 1; a(0,1) = 2; a(1,0) = 3; a(1,1) = 4;
    bool test_condition = (deter(a) == -2);
    if ( test_condition ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }

  ////////////////////////////////////////////////////////////////////////////////
  void check_determinant_3D()
  {
    Matrix<double, 3, 3> a;
    a(0,0) = 1; a(0,1) = 2;  a(0,2) = 3;
    a(1,0) = 4; a(1,1) = 15; a(1,2) = 6;
    a(2,0) = 7; a(2,1) = 8;  a(2,2) = 9;
    bool test_condition = (deter(a) == -120);
    if ( test_condition ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }

  ////////////////////////////////////////////////////////////////////////////////
  void check_con()
  {
    Matrix<double, 3, 3> a;
    a(0,0) = 1; a(0,1) = 2;  a(0,2) = 3;
    a(1,0) = 4; a(1,1) = 15; a(1,2) = 6;
    a(2,0) = 7; a(2,1) = 8;  a(2,2) = 9;
    Matrix<double, 3, 3> b = 5; b += a;
    bool test_condition = (con(a,b) == 760);
    if ( test_condition ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }

  ////////////////////////////////////////////////////////////////////////////////
  void check_con2()
  {
    Matrix<double, 3, 3> a;
    a(0,0) = 1; a(0,1) = 2;  a(0,2) = 3;
    a(1,0) = 4; a(1,1) = 15; a(1,2) = 6;
    a(2,0) = 7; a(2,1) = 8;  a(2,2) = 9;
    Matrix<double, 3, 3> b = 5; b += a;
    bool test_condition = (con2(a,b) == 760);
    if ( test_condition ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }

  ////////////////////////////////////////////////////////////////////////////////
  void check_rowp1()
  {
    Matrix<double, 3, 3> a;
    a(0,0) = 1; a(0,1) = 2;  a(0,2) = 3;
    a(1,0) = 4; a(1,1) = 15; a(1,2) = 6;
    a(2,0) = 7; a(2,1) = 8;  a(2,2) = 9;
    Vector<double,2> b = rowp1( 0, a );
    bool test_condition = (b(0) == 2 && b(1) == 3);
    b = rowp1( 1, a );
    test_condition = test_condition && (b(0) == 15 && b(1) == 6);
    b = rowp1( 2, a );
    test_condition = test_condition && (b(0) == 8 && b(1) == 9);
    if ( test_condition ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }

  ////////////////////////////////////////////////////////////////////////////////
  void check_colp1()
  {
    Matrix<double, 3, 3> a;
    a(0,0) = 1; a(0,1) = 2;  a(0,2) = 3;
    a(1,0) = 4; a(1,1) = 15; a(1,2) = 6;
    a(2,0) = 7; a(2,1) = 8;  a(2,2) = 9;
    Vector<double,2> b = colp1( 0, a );
    bool test_condition = (b(0) == 4 && b(1) == 7);
    b = colp1( 1, a );
    test_condition = test_condition && (b(0) == 15 && b(1) == 8);
    b = colp1( 2, a );
    test_condition = test_condition && (b(0) == 6 && b(1) == 9);
    if ( test_condition ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }

  ////////////////////////////////////////////////////////////////////////////////
  void check_mini()
  {
    Matrix<double, 3, 3> a;
    a(0,0) = 1; a(0,1) = 2;  a(0,2) = 3;
    a(1,0) = 4; a(1,1) = 15; a(1,2) = 6;
    a(2,0) = 7; a(2,1) = 8;  a(2,2) = 9;
    Matrix<double, 2, 2> b;
    b(0,0) = -10; b(0,1) = -11;
    b(1,0) = -12; b(1,1) = -13;
    mini( b, a );
    bool test_condition = (b(0,0) == 15 && b(0,1) == 6
                            && b(1,0) == 8 && b(1,1) == 9);
    if ( test_condition ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }

  ////////////////////////////////////////////////////////////////////////////////
  void check_tmini()
  {
    Matrix<double, 3, 3> a;
    a(0,0) = 1; a(0,1) = 2;  a(0,2) = 3;
    a(1,0) = 4; a(1,1) = 15; a(1,2) = 6;
    a(2,0) = 7; a(2,1) = 8;  a(2,2) = 9;
    Matrix<double, 2, 2> b;
    b(0,0) = -10; b(0,1) = -11;
    b(1,0) = -12; b(1,1) = -13;
    tmini( a, b );
    bool test_condition = (a(0,0) == 1 && a(0,1) == 2 && a(0,2) == 3
                            && a(1,0) == 4 && a(1,1) == -10 && a(1,2) == -11
                            && a(2,0) == 7 && a(2,1) == -12 && a(2,2) == -13);
    if ( test_condition ) passed( __FUNCTION__ ); else failed( __FUNCTION__ );
  }


}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void run_matrix_tests()
{
  cout << "Running matrix checks" << endl;
  using namespace matrix_checks;

  check_copy();
  check_assignment();
  check_inplace_addition();
  check_inplace_subtraction();
  check_inplace_scalar_multiplication();
  check_negation();
  check_addition();
  check_subtraction();
  check_scalar_multiplication();

  check_matrix_matrix_multiplication();
  check_vector_matrix_multiplication();
  check_matrix_vector_multiplication();

  check_transpose();
  check_determinant_2D();
  check_determinant_3D();
  check_con();
  check_con2();

  check_rowp1();
  check_colp1();
  check_mini();
  check_tmini();

  cout << endl;
}

