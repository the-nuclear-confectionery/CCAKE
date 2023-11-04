#include <gtest/gtest.h>
#include "kernel.h"
#include "vector.h"

#define TOL 1.E-14

TEST( KernelTest, Kernel1D_q_0p5 )
{
  double dist = .5;
  double h = 1.;
  EXPECT_LT( fabs(ccake::SPHkernel<1>::kernel(dist,h)-23./48.) , TOL )
  << "Kernel value is "<< ccake::SPHkernel<1>::kernel(dist,h) 
  << " instead of expected "  << 23./48.;
}

TEST( KernelTest, Kernel1D_q_p10 )
{
  double dist = 1.;
  double h = 1.;
  EXPECT_LT( fabs(ccake::SPHkernel<1>::kernel(dist,h)-1./6.) , TOL )
  << "Kernel value is "<< ccake::SPHkernel<1>::kernel(dist,h) 
  << " instead of expected "  << 1./6.;
}

TEST( KernelTest, Kernel1D_q_1p5 )
{
  double dist = 1.5;
  double h = 1;
  EXPECT_LT( fabs(ccake::SPHkernel<1>::kernel(dist,h)-1./48.) , TOL )
  << "Kernel value is "<< ccake::SPHkernel<1>::kernel(dist,h) 
  << " instead of expected "  << 1./48.;
}

TEST( KernelTest, Kernel1D_q_2p0 )
{
  double dist = 2.0;
  double h = 1;
  EXPECT_LT( fabs(ccake::SPHkernel<1>::kernel(dist,h)) , TOL )
  << "Kernel value is "<< ccake::SPHkernel<1>::kernel(dist,h) 
  << " instead of expected "  << 0.;
}


TEST( KernelTest, Kernel2D_q_0p5 )
{
  double dist = .5;
  double h = 1.;
  EXPECT_LT( fabs(ccake::SPHkernel<2>::kernel(dist,h)-115./112./M_PI) , TOL )
  << "Kernel value is "<< ccake::SPHkernel<2>::kernel(dist,h) 
  << " instead of expected "  << 115./112./M_PI;
}

TEST( KernelTest, Kernel2D_q_1p0 )
{
  double dist = 1.;
  double h = 1.;
  EXPECT_LT( fabs(ccake::SPHkernel<2>::kernel(dist,h)-5./14./M_PI) , TOL )
  << "Kernel value is "<< ccake::SPHkernel<2>::kernel(dist,h) 
  << " instead of expected "  << 115./112./M_PI;
}

TEST( KernelTest, Kernel2D_q_1p5 )
{
  double dist = 1.5;
  double h = 1;
  EXPECT_LT( fabs(ccake::SPHkernel<2>::kernel(dist,h)-5./112./M_PI) , TOL )
  << "Kernel value is "<< ccake::SPHkernel<2>::kernel(dist,h) 
  << " instead of expected "  << 5./112./M_PI;
}

TEST( KernelTest, Kernel2D_q_2p0 )
{
  double dist = 2.0;
  double h = 1;
  EXPECT_LT( fabs(ccake::SPHkernel<2>::kernel(dist,h)) , TOL )
  << "Kernel value is "<< ccake::SPHkernel<2>::kernel(dist,h) 
  << " instead of expected "  << 0.;
}


TEST( KernelTest, Kernel3D_q_0p5 )
{
  double dist = .5;
  double h = 1.;
  EXPECT_LT( fabs(ccake::SPHkernel<3>::kernel(dist,h)-23./32./M_PI) , TOL )
  << "Kernel value is "<< ccake::SPHkernel<3>::kernel(dist,h) 
  << " instead of expected "  << 23./32./M_PI;
}

TEST( KernelTest, Kernel3D_q_1p0 )
{
  double dist = 1.;
  double h = 1.;
  EXPECT_LT( fabs(ccake::SPHkernel<3>::kernel(dist,h)-1./4./M_PI) , TOL )
  << "Kernel value is "<< ccake::SPHkernel<3>::kernel(dist,h) 
  << " instead of expected "  << 1./4./M_PI;
}

TEST( KernelTest, Kernel3D_q_1p5 )
{
  double dist = 1.5;
  double h = 1;
  EXPECT_LT( fabs(ccake::SPHkernel<3>::kernel(dist,h)-1./32./M_PI) , TOL )
  << "Kernel value is "<< ccake::SPHkernel<3>::kernel(dist,h) 
  << " instead of expected "  << 1./32./M_PI;
}

TEST( KernelTest, Kernel3D_q_2p0 )
{
  double dist = 2.0;
  double h = 1;
  EXPECT_LT( fabs(ccake::SPHkernel<3>::kernel(dist,h)) , TOL )
  << "Kernel value is "<< ccake::SPHkernel<3>::kernel(dist,h) 
  << " instead of expected "  << 0.;
}