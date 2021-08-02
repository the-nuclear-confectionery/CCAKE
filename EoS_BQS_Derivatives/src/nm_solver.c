#include <math.h>
#include <stdio.h>
#define NRANSI

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

#include "../include/Functions.h"
#include "../include/Variables.h"
#include "../include/nm_solver.h"

void solve ( double densities[], double sols[] )
{
	double e0 = densities[0], B0 = densities[1], S0 = densities[2], Q0 = densities[3];
	double Tout = sols[0], muBout = sols[1], muSout = sols[2], muQout = sols[3];

	const int maxTries = 10;
	const double ACCURACY = 1e-6;
	const double hbarc = 197.327;
	const double hbarc3 = hbarc*hbarc*hbarc;

	//initial guess
	//Tout = hbarc; muBout = 0.0; muSout = 0.0; muQout = 0.0;
	double T2 = Tout*Tout;
	double T3 = T2*Tout;
	double T4 = T3*Tout;

	// compute initial estimates
	double Plocal = T4*PressTaylor(Tout, muBout, muQout, muSout)/hbarc3;
	double slocal = T3*EntrTaylor(Tout, muBout, muQout, muSout)/hbarc3;
	double Bsol = T3*BarDensTaylor(Tout, muBout, muQout, muSout)/hbarc3;
	double Ssol = T3*StrDensTaylor(Tout, muBout, muQout, muSout)/hbarc3;
	double Qsol = T3*ChDensTaylor(Tout, muBout, muQout, muSout)/hbarc3;
	double esol = slocal*Tout - Plocal + muBout*Bsol + muQout*Qsol + muSout*Ssol;

	printf("e0 = %15.8f\n", e0);
	printf("B0 = %15.8f\n", B0);
	printf("S0 = %15.8f\n", S0);
	printf("Q0 = %15.8f\n", Q0);
	printf("esol = %15.8f\n", esol);
	printf("Bsol = %15.8f\n", Bsol);
	printf("Ssol = %15.8f\n", Ssol);
	printf("Qsol = %15.8f\n\n", Qsol);

	//bool not_converged = abs(esol-e0) > ACCURACY or abs(Bsol-B0) > ACCURACY
	//					  or abs(Ssol-S0) > ACCURACY or abs(Qsol-Q0) > ACCURACY;

	int iter = 0;
	while ( (fabs(esol-e0) > ACCURACY || fabs(Bsol-B0) > ACCURACY
			  || fabs(Ssol-S0) > ACCURACY || fabs(Qsol-Q0) > ACCURACY) && iter++ < maxTries )
	{
		printf("iter = %5d\n", iter);
		T2 = Tout*Tout; T3 = T2*Tout;

		double dBdT   = T2*P2TB(Tout, muBout, muQout, muSout);
		double dBdmuB = T2*P2B2(Tout, muBout, muQout, muSout);
		double dBdmuS = T2*P2BS(Tout, muBout, muQout, muSout);
		double dBdmuQ = T2*P2BQ(Tout, muBout, muQout, muSout);

		double dSdT   = T2*P2TS(Tout, muBout, muQout, muSout);
		double dSdmuB = T2*P2BS(Tout, muBout, muQout, muSout);
		double dSdmuS = T2*P2S2(Tout, muBout, muQout, muSout);
		double dSdmuQ = T2*P2QS(Tout, muBout, muQout, muSout);

		double dQdT   = T2*P2TQ(Tout, muBout, muQout, muSout);
		double dQdmuB = T2*P2BQ(Tout, muBout, muQout, muSout);
		double dQdmuS = T2*P2QS(Tout, muBout, muQout, muSout);
		double dQdmuQ = T2*P2Q2(Tout, muBout, muQout, muSout);

		double dedT   = T3*P2T2(Tout, muBout, muQout, muSout)
						+ muBout*dBdT   + muSout*dSdT   + muQout*dQdT;
		double dedmuB = T3*P2TB(Tout, muBout, muQout, muSout)
						+ muBout*dBdmuB + muSout*dSdmuB + muQout*dQdmuB;
		double dedmuS = T3*P2TS(Tout, muBout, muQout, muSout)
						+ muBout*dBdmuS + muSout*dSdmuS + muQout*dQdmuS;
		double dedmuQ = T3*P2TQ(Tout, muBout, muQout, muSout)
						+ muBout*dBdmuQ + muSout*dSdmuQ + muQout*dQdmuQ;

		double a_data[] = { dedT/hbarc3, dedmuB/hbarc3, dedmuS/hbarc3, dedmuQ/hbarc3,
							dBdT/hbarc3, dBdmuB/hbarc3, dBdmuS/hbarc3, dBdmuQ/hbarc3,
							dSdT/hbarc3, dSdmuB/hbarc3, dSdmuS/hbarc3, dSdmuQ/hbarc3,
							dQdT/hbarc3, dQdmuB/hbarc3, dQdmuS/hbarc3, dQdmuQ/hbarc3 };
		double b_data[] = { esol-e0, Bsol-B0, Ssol-S0, Qsol-Q0 };
		
		gsl_matrix_view m = gsl_matrix_view_array (a_data, 4, 4);
		gsl_vector_view b = gsl_vector_view_array (b_data, 4);
		gsl_vector *x = gsl_vector_alloc (4);
		
		int s;
		gsl_permutation * p = gsl_permutation_alloc (4);
		gsl_linalg_LU_decomp (&m.matrix, p, &s);
		gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
		
//		esol -= gsl_vector_get(x, 0);
//		Bsol -= gsl_vector_get(x, 1);
//		Ssol -= gsl_vector_get(x, 2);
//		Qsol -= gsl_vector_get(x, 3);
		Tout -= gsl_vector_get(x, 0);
		muBout -= gsl_vector_get(x, 1);
		muSout -= gsl_vector_get(x, 2);
		muQout -= gsl_vector_get(x, 3);

	printf("e0 = %15.8f\n", e0);
	printf("B0 = %15.8f\n", B0);
	printf("S0 = %15.8f\n", S0);
	printf("Q0 = %15.8f\n", Q0);
	printf("esol = %15.8f\n", esol);
	printf("Bsol = %15.8f\n", Bsol);
	printf("Ssol = %15.8f\n", Ssol);
	printf("Qsol = %15.8f\n", Qsol);

		printf("%15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f"
				"%15.12f %15.12f %15.12f %15.12f %15.12f\n",
				esol, e0, fabs(esol-e0), Bsol, B0, fabs(Bsol-B0),
				Ssol, S0, fabs(Ssol-S0), Qsol, Q0, fabs(Qsol-Q0), ACCURACY);

		
	printf("********************************************************************************\n\n");

		gsl_permutation_free (p);
		gsl_vector_free (x);

		//not_converged = abs(esol-e0) > ACCURACY or abs(Bsol-B0) > ACCURACY
		//				or abs(Ssol-S0) > ACCURACY or abs(Qsol-Q0) > ACCURACY;

	}

	// otherwise, store some other combination that's maybe nearby
	if ( fabs(esol-e0) > ACCURACY || fabs(Bsol-B0) > ACCURACY
		 || fabs(Ssol-S0) > ACCURACY || fabs(Qsol-Q0) > ACCURACY )
	{
		printf("Did not find a solution!\n");
		printf("%15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f"
				"%15.12f %15.12f %15.12f %15.12f %15.12f\n",
				esol, e0, fabs(esol-e0), Bsol, B0, fabs(Bsol-B0),
				Ssol, S0, fabs(Ssol-S0), Qsol, Q0, fabs(Qsol-Q0), ACCURACY);

		e0 = esol;
		B0 = Bsol;
		S0 = Ssol;
		Q0 = Qsol;
	}

	densities[0] = e0; densities[1] = B0; densities[2] = S0; densities[3] = Q0;
	sols[0] = Tout; sols[1] = muBout; sols[2] = muSout; sols[3] = muQout;

	return;
}

#undef NRANSI

