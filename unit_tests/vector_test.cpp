#include <gtest/gtest.h>
#include "vector.h"

#define TOL 1.E-14

TEST( VectorTest, Assignment1D )
{
  Vector<int,1> a = {1};
  Vector<int,1> b = a;
  EXPECT_EQ( a(0), b(0) );
}

TEST( VectorTest, Assignment2D )
{
  Vector<int,2> a2 = {2,3};
  Vector<int,2> b2 = a2;
  EXPECT_EQ( a2(0), b2(0) );
  EXPECT_EQ( a2(1), b2(1) );
}

TEST( VectorTest, Assignment3D )
{
  Vector<int,3> a3 = {3,4,5};
  Vector<int,3> b3 = a3;
  EXPECT_EQ( a3(0), b3(0) );
  EXPECT_EQ( a3(1), b3(1) );
  EXPECT_EQ( a3(2), b3(2) );
}

TEST ( VectorTest, AssignmentScalar )
{
  Vector<int,3> a = {1,2,3};
  double b = 7;
  a = b;
  for (int i=0; i<3; i++)
    EXPECT_EQ( a(i), 7 );
}

TEST( VectorTest, AssignmentAdd)
{
  Vector<int,3> a = {1,2,3};
  Vector<int,3> a_original = {1,2,3};
  Vector<int,3> b = {4,5,6};
  a += b;
  for (int i=0; i<3; i++)
    EXPECT_EQ( a(i), a_original(i) + b(i) );
}

TEST( VectorTest, AssignmentSubtract)
{
  Vector<int,3> a = {1,2,3};
  Vector<int,3> a_original = {1,2,3};
  Vector<int,3> b = {4,5,6};
  a -= b;
  for (int i=0; i<3; i++)
    EXPECT_EQ( a(i), a_original(i) - b(i) );
}

TEST( VectorTest, AssignmentMultiply)
{
  Vector<int,3> a = {1,2,3};
  Vector<int,3> a_original = {1,2,3};
  double b = 7;
  a *= b;
  for (int i=0; i<3; i++)
    EXPECT_EQ( a(i), a_original(i) * b );
}

TEST ( VectorTest, Addition )
{
  Vector<int,3> a = {1,2,3};
  Vector<int,3> b = {4,5,6};
  Vector<int,3> c = a + b;
  for (int i=0; i<3; i++)
    EXPECT_EQ( c(i), a(i) + b(i) );
}


TEST ( VectorTest, SignReverseInt )
{
  Vector<int,3> a = {1,2,3};
  Vector<int,3> b = -a;
  for (int i=0; i<3; i++)
    EXPECT_EQ( b(i), -a(i) );
}

TEST ( VectorTest, SignReverseDouble )
{
  Vector<double,3> a = {1,2,3};
  Vector<double,3> b = -a;
  for (int i=0; i<3; i++)
    EXPECT_LT( fabs(b(i)+a(i)), TOL );
}

TEST ( VectorTest, Subtraction )
{
  Vector<int,3> a = {1,2,3};
  Vector<int,3> b = {4,5,6};
  Vector<int,3> c = a - b;
  for (int i=0; i<3; i++)
    EXPECT_EQ( c(i), a(i) - b(i) );
}

TEST ( VectorTest, ScalarMultiplication )
{
  Vector<int,3> a = {1,2,3};
  int b = 7;
  Vector<int,3> c = b * a;
  for (int i=0; i<3; i++)
    EXPECT_EQ( c(i), a(i) * b );
}

TEST( VectorTest, InnerProduct )
{
  Vector<int,3> a = {1,2,3};
  Vector<int,3> b = {4,5,6};
  int c = inner(a,b);
  EXPECT_EQ( c, 32 );
}

TEST( VectorTest, Norm )
{
  Vector<int,2> a = {3,4};
  double b = Norm(a);
  EXPECT_LT( fabs(b - 5), TOL );
}

TEST( VectorTest, Norm2 )
{
  Vector<int,2> a = {3,4};
  double b = Norm2(a);
  EXPECT_LT( fabs(b - 25), TOL );
}
