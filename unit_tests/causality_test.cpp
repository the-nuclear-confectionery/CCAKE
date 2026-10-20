#include <gtest/gtest.h>
#include <cmath>
#include <random>

#include "milne.hpp"          // eom_causality.hpp is templated on the Matrix type
#include "eom_causality.hpp"

/// @file causality_test.cpp
/// @brief Regression tests for the shear-eigenvalue extraction behind the BDNR
/// causality check (eom_causality.hpp).

using namespace ccake::causality;

//==============================================================================
// Lambda extraction
//==============================================================================

/// pi^mu_nu is traceless, so the extracted eigenvalues must sum to zero -- by
/// construction of the depressed cubic, not merely to round-off luck.
TEST( CausalityLambda, Traceless )
{
  std::mt19937_64 rng(20260716);
  std::uniform_real_distribution<double> u(-0.5, 0.5);

  for ( int trial = 0; trial < 200; ++trial )
  {
    const double t2 = 1.0 + 3.0*std::fabs(u(rng));

    // a symmetric pi^{mu nu} that is traceless w.r.t. g = diag(1,-1,-1,-t2):
    // pi^33 = (pi^00 - pi^11 - pi^22)/t2, exactly as reset_pi_tensor builds it.
    ccake::milne::Matrix<double,4,4> shv;
    for ( int i = 0; i < 4; ++i )
      for ( int j = 0; j <= i; ++j )
      {
        const double val = u(rng);
        shv(i,j) = val; shv(j,i) = val;
      }
    shv(3,3) = (shv(0,0) - shv(1,1) - shv(2,2))/t2;

    double Lam[3];
    shear_eigenvalues(shv, t2, Lam);

    EXPECT_NEAR( Lam[0] + Lam[1] + Lam[2], 0.0, 1.0e-10 ) << "trial " << trial;
    EXPECT_LE( Lam[0], Lam[1] );
    EXPECT_LE( Lam[1], Lam[2] );
  }
}

/// Pins the metric-signature convention, which is the subtlest part of this code.
///
/// The paper is mostly-plus, CCAKE is mostly-minus. The contravariant pi^{ij} is the
/// same physical object in both (it is the deviation of T^{ij} from the isotropic
/// part), but the MIXED tensor pi^mu_nu flips sign between them -- so the paper's
/// Lambda_a are MINUS the eigenvalues of CCAKE's mixed tensor. Equivalently, and
/// this is what we assert: in the local rest frame the paper's Lambda_a are just the
/// eigenvalues of the physical spatial shear stress pi^{ij}.
///
/// Get this backwards and every condition, old and new, is silently evaluated on
/// -Lambda. Nothing else in the suite would catch it, because the conditions are not
/// symmetric under Lambda -> -Lambda but they are all still *computable*.
TEST( CausalityLambda, RestFrameSignConvention )
{
  const double t2 = 1.0;   // tau = 1, so the eta-eta slot carries no extra factor

  ccake::milne::Matrix<double,4,4> shv;
  for ( int i = 0; i < 4; ++i ) for ( int j = 0; j < 4; ++j ) shv(i,j) = 0.0;

  // A diagonal rest-frame shear stress. Traceless w.r.t. g = diag(1,-1,-1,-1):
  //   g_{mn} pi^{mn} = 0 - 0.3 - 0.5 + 0.8 = 0
  shv(1,1) =  0.3;
  shv(2,2) =  0.5;
  shv(3,3) = -0.8;

  double Lam[3];
  shear_eigenvalues(shv, t2, Lam);

  // the eigenvalues of pi^{ij} itself, ascending -- NOT their negatives
  EXPECT_NEAR( Lam[0], -0.8, 1.0e-12 );
  EXPECT_NEAR( Lam[1],  0.3, 1.0e-12 );
  EXPECT_NEAR( Lam[2],  0.5, 1.0e-12 );
}

/// Vanishing shear stress must give Lambda = 0 rather than a division blow-up.
TEST( CausalityLambda, ZeroShear )
{
  ccake::milne::Matrix<double,4,4> shv;
  for ( int i = 0; i < 4; ++i ) for ( int j = 0; j < 4; ++j ) shv(i,j) = 0.0;

  double Lam[3];
  shear_eigenvalues(shv, 1.0, Lam);
  for ( int k = 0; k < 3; ++k ) EXPECT_DOUBLE_EQ( Lam[k], 0.0 );
}
