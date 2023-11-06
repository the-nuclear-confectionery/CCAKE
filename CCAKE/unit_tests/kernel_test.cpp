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

TEST( GradKernelTest, GradKernel2D_q_0p5 )
{
  double pos_a[2] = {.25,.0};
  double pos_b[2] = {-.25,.0};
  double rel_pos[2];
  double const target_val = -75./56./M_PI;
  double h = 1;
  double gradK[2];
  
  double dist = ccake::SPHkernel<2>::distance(pos_a,pos_b);
  for (int idir=0; idir<2;++idir) rel_pos[idir] = pos_a[idir] - pos_b[idir];
  ccake::SPHkernel<2>::gradKernel(rel_pos,dist,h,gradK);
  
  EXPECT_LT( fabs(gradK[0] - target_val), TOL) << "gradK[0] = " << gradK[0] <<". Expected " << target_val;
  EXPECT_LT( fabs(gradK[1]), TOL) << "gradK[1] = " << gradK[1] <<". Expected 0" << endl;
  
  double pos_a2[2] = {0.,.25};
  double pos_b2[2] = {0.,-.25};
  for (int idir=0; idir<2;++idir) rel_pos[idir] = pos_a2[idir] - pos_b2[idir];
  dist = ccake::SPHkernel<2>::distance(pos_a2,pos_b2);
  ccake::SPHkernel<2>::gradKernel(rel_pos,dist,h,gradK);
  EXPECT_LT( fabs(gradK[0]), TOL) << "gradK[0] = " << gradK[0] <<". Expected 0" << endl;
  EXPECT_LT( fabs(gradK[1] - target_val), TOL) << "gradK[1] = " << gradK[1] <<". Expected " << target_val;
  
}

TEST( GradKernelTest, GradKernel2D_q_1p0 )
{
  double pos_a[2] = {.5,.0};
  double pos_b[2] = {-.5,.0};
  double rel_pos[2];
  double const target_val = -15./14./M_PI;
  double h = 1;
  double gradK[2];
  
  double dist = ccake::SPHkernel<2>::distance(pos_a,pos_b);
  for (int idir=0; idir<2;++idir) rel_pos[idir] = pos_a[idir] - pos_b[idir];
  ccake::SPHkernel<2>::gradKernel(rel_pos,dist,h,gradK);
  
  EXPECT_LT( fabs(gradK[0] - target_val), TOL) << "gradK[0] = " << gradK[0] <<". Expected " << target_val;
  EXPECT_LT( fabs(gradK[1]), TOL) << "gradK[1] = " << gradK[1] <<". Expected 0" << endl;
  
  double pos_a2[2] = {0.,.5};
  double pos_b2[2] = {0.,-.5};
  for (int idir=0; idir<2;++idir) rel_pos[idir] = pos_a2[idir] - pos_b2[idir];
  dist = ccake::SPHkernel<2>::distance(pos_a2,pos_b2);
  ccake::SPHkernel<2>::gradKernel(rel_pos,dist,h,gradK);
  EXPECT_LT( fabs(gradK[0]), TOL) << "gradK[0] = " << gradK[0] <<". Expected 0" << endl;
  EXPECT_LT( fabs(gradK[1] - target_val), TOL) << "gradK[1] = " << gradK[1] <<". Expected " << target_val;
  
}

TEST( GradKernelTest, GradKernel2D_q_1p5 )
{
  double pos_a[2] = {.75,.0};
  double pos_b[2] = {-.75,.0};
  double rel_pos[2];
  double const target_val = -15./56./M_PI;
  double h = 1;
  double gradK[2];
  
  double dist = ccake::SPHkernel<2>::distance(pos_a,pos_b);
  for (int idir=0; idir<2;++idir) rel_pos[idir] = pos_a[idir] - pos_b[idir];
  ccake::SPHkernel<2>::gradKernel(rel_pos,dist,h,gradK);
  
  EXPECT_LT( fabs(gradK[0] - target_val), TOL) << "gradK[0] = " << gradK[0] <<". Expected " << target_val;
  EXPECT_LT( fabs(gradK[1]), TOL) << "gradK[1] = " << gradK[1] <<". Expected 0" << endl;
  
  double pos_a2[2] = {0.,.75};
  double pos_b2[2] = {0.,-.75};
  for (int idir=0; idir<2;++idir) rel_pos[idir] = pos_a2[idir] - pos_b2[idir];
  dist = ccake::SPHkernel<2>::distance(pos_a2,pos_b2);
  ccake::SPHkernel<2>::gradKernel(rel_pos,dist,h,gradK);
  EXPECT_LT( fabs(gradK[0]), TOL) << "gradK[0] = " << gradK[0] <<". Expected 0" << endl;
  EXPECT_LT( fabs(gradK[1] - target_val), TOL) << "gradK[1] = " << gradK[1] <<". Expected " << target_val;
  
}