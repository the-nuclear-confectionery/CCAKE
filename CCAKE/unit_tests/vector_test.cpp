#include <gtest/gtest.h>
#include "vector.h"

TEST( VectorTest, Assignment1D )
{
  Vector<int,3> a = {1};
  Vector<int,3> b = a;
  EXPECT_EQ( a(0), b(0) );
}

TEST( VectorTest, Assignment2D )
{
  Vector<int,3> a2 = {2,3};
  Vector<int,3> b2 = a2;
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