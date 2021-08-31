#include <math.h>
#include <stdio.h>
#define NRANSI

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

#include "../include/Functions.h"
#include "../include/Variables.h"
#include "../include/nm_solver.h"

double Tmin = 30.0, Tmax = 800.0;
double amuBmax = 450.0, amuSmax = 450.0, amuQmax = 450.0;

void solve ( double densities[], double sols[] )
{
	double eTarget = densities[0], BTarget = densities[1], STarget = densities[2], QTarget = densities[3];
	double Tout = sols[0], muBout = sols[1], muSout = sols[2], muQout = sols[3];

	const int maxTries = 10000;
	const double ACCURACY = 1e-4;
	const double hbarc = 197.327;
	const double hbarc3 = hbarc*hbarc*hbarc;

	//initial guess
	//Tout = hbarc; muBout = 0.0; muSout = 0.0; muQout = 0.0;
	double T2 = Tout*Tout;
	double T3 = T2*Tout;
	double T4 = T3*Tout;

	// compute initial estimates
	double Plocal = T4*PressTaylor(Tout, muBout, muQout, muSout)/hbarc3;			// MeV/fm^3
	double slocal = T3*EntrTaylor(Tout, muBout, muQout, muSout)/hbarc3;				// 1/fm^3
	double Bsol = T3*BarDensTaylor(Tout, muBout, muQout, muSout)/hbarc3;			// 1/fm^3
	double Ssol = T3*StrDensTaylor(Tout, muBout, muQout, muSout)/hbarc3;			// 1/fm^3
	double Qsol = T3*ChDensTaylor(Tout, muBout, muQout, muSout)/hbarc3;				// 1/fm^3
	double esol = slocal*Tout - Plocal + muBout*Bsol + muQout*Qsol + muSout*Ssol;	// MeV/fm^3

	printf("eTarget = %15.8f\n", eTarget);
	printf("BTarget = %15.8f\n", BTarget);
	printf("STarget = %15.8f\n", STarget);
	printf("QTarget = %15.8f\n", QTarget);
	printf("esol = %15.8f\n", esol);
	printf("Bsol = %15.8f\n", Bsol);
	printf("Ssol = %15.8f\n", Ssol);
	printf("Qsol = %15.8f\n\n", Qsol);

	int iter = 0;
	while ( (fabs(esol-eTarget) > ACCURACY*eTarget
			|| fabs(Bsol-BTarget) > fmax(ACCURACY*fabs(BTarget), ACCURACY)
			|| fabs(Ssol-STarget) > fmax(ACCURACY*fabs(STarget), ACCURACY)
			|| fabs(Qsol-QTarget) > fmax(ACCURACY*fabs(QTarget), ACCURACY))
			&& iter++ < maxTries )
	{
		printf("iter = %5d\n", iter);
		T2 = Tout*Tout; T3 = T2*Tout; T4 = T3*Tout;

		double dBdT   = T2*P2TB(Tout, muBout, muQout, muSout);				// MeV^2
		double dBdmuB = T2*P2B2(Tout, muBout, muQout, muSout);				// MeV^2
		double dBdmuS = T2*P2BS(Tout, muBout, muQout, muSout);				// MeV^2
		double dBdmuQ = T2*P2BQ(Tout, muBout, muQout, muSout);				// MeV^2

		double dSdT   = T2*P2TS(Tout, muBout, muQout, muSout);				// MeV^2
		double dSdmuB = T2*P2BS(Tout, muBout, muQout, muSout);				// MeV^2
		double dSdmuS = T2*P2S2(Tout, muBout, muQout, muSout);				// MeV^2
		double dSdmuQ = T2*P2QS(Tout, muBout, muQout, muSout);				// MeV^2

		double dQdT   = T2*P2TQ(Tout, muBout, muQout, muSout);				// MeV^2
		double dQdmuB = T2*P2BQ(Tout, muBout, muQout, muSout);				// MeV^2
		double dQdmuS = T2*P2QS(Tout, muBout, muQout, muSout);				// MeV^2
		double dQdmuQ = T2*P2Q2(Tout, muBout, muQout, muSout);				// MeV^2

		double dedT   = T3*P2T2(Tout, muBout, muQout, muSout)
						+ muBout*dBdT   + muSout*dSdT   + muQout*dQdT;		// MeV^3
		double dedmuB = T3*P2TB(Tout, muBout, muQout, muSout)
						+ muBout*dBdmuB + muSout*dSdmuB + muQout*dQdmuB;	// MeV^3
		double dedmuS = T3*P2TS(Tout, muBout, muQout, muSout)
						+ muBout*dBdmuS + muSout*dSdmuS + muQout*dQdmuS;	// MeV^3
		double dedmuQ = T3*P2TQ(Tout, muBout, muQout, muSout)
						+ muBout*dBdmuQ + muSout*dSdmuQ + muQout*dQdmuQ;	// MeV^3

		double a_data[] = { dedT/hbarc3, dedmuB/hbarc3, dedmuS/hbarc3, dedmuQ/hbarc3,
							dBdT/hbarc3, dBdmuB/hbarc3, dBdmuS/hbarc3, dBdmuQ/hbarc3,
							dSdT/hbarc3, dSdmuB/hbarc3, dSdmuS/hbarc3, dSdmuQ/hbarc3,
							dQdT/hbarc3, dQdmuB/hbarc3, dQdmuS/hbarc3, dQdmuQ/hbarc3 };
		double b_data[] = { esol-eTarget, Bsol-BTarget, Ssol-STarget, Qsol-QTarget };

		//for (int ii = 0; ii < 16; ii++) printf("a_data[%5d] = %15.12f\n", ii, a_data[ii]);
		//for (int ii = 0; ii < 4; ii++) printf("b_data[%5d] = %15.12f\n", ii, b_data[ii]);

		gsl_matrix_view m = gsl_matrix_view_array (a_data, 4, 4);
		gsl_vector_view b = gsl_vector_view_array (b_data, 4);
		gsl_vector *x = gsl_vector_alloc (4);
		
		int s;
		gsl_permutation * p = gsl_permutation_alloc (4);
		gsl_linalg_LU_decomp (&m.matrix, p, &s);
		gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
		
		// if out of range
		if ( Tout - gsl_vector_get(x, 0) < Tmin
				|| Tout - gsl_vector_get(x, 0) > Tmax
				|| abs(muBout - gsl_vector_get(x, 1)) > amuBmax
				|| abs(muSout - gsl_vector_get(x, 2)) > amuSmax
				|| abs(muQout - gsl_vector_get(x, 3)) > amuQmax )
		{
			sols[0] = -1.0;	// indicates failure since T >= 0
			return;			// exit prematurely
		}

		// otherwise, continue
		Tout -= gsl_vector_get(x, 0);
		muBout -= gsl_vector_get(x, 1);
		muSout -= gsl_vector_get(x, 2);
		muQout -= gsl_vector_get(x, 3);

		printf("Tout = %15.8f\n", Tout);
		printf("muBout = %15.8f\n", muBout);
		printf("muSout = %15.8f\n", muSout);
		printf("muQout = %15.8f\n", muQout);

		// update previous estimates
		Plocal = T4*PressTaylor(Tout, muBout, muQout, muSout)/hbarc3;
		slocal = T3*EntrTaylor(Tout, muBout, muQout, muSout)/hbarc3;
		Bsol = T3*BarDensTaylor(Tout, muBout, muQout, muSout)/hbarc3;
		Ssol = T3*StrDensTaylor(Tout, muBout, muQout, muSout)/hbarc3;
		Qsol = T3*ChDensTaylor(Tout, muBout, muQout, muSout)/hbarc3;
		esol = slocal*Tout - Plocal + muBout*Bsol + muQout*Qsol + muSout*Ssol;

		printf("Plocal = %15.8f\n", Plocal);
		printf("slocal = %15.8f\n", slocal);
	
		printf("eTarget = %15.8f\n", eTarget);
		printf("BTarget = %15.8f\n", BTarget);
		printf("STarget = %15.8f\n", STarget);
		printf("QTarget = %15.8f\n", QTarget);
		printf("esol = %15.8f\n", esol);
		printf("Bsol = %15.8f\n", Bsol);
		printf("Ssol = %15.8f\n", Ssol);
		printf("Qsol = %15.8f\n", Qsol);
	
		printf("%15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f"
				"%15.12f %15.12f %15.12f %15.12f %15.12f\n",
				esol, eTarget, fabs(esol-eTarget), Bsol, BTarget, fabs(Bsol-BTarget),
				Ssol, STarget, fabs(Ssol-STarget), Qsol, QTarget, fabs(Qsol-QTarget), ACCURACY);
	
		printf("********************************************************************************\n\n");

		gsl_permutation_free (p);
		gsl_vector_free (x);

	}

	if ( fabs(esol-eTarget) > ACCURACY || fabs(Bsol-BTarget) > ACCURACY
		 || fabs(Ssol-STarget) > ACCURACY || fabs(Qsol-QTarget) > ACCURACY )
		{
			sols[0] = -1.0;	// indicates failure since T >= 0
			return;			// exit prematurely
		}

	printf("Final iter = %d\n", iter);

	sols[0] = Tout; sols[1] = muBout; sols[2] = muSout; sols[3] = muQout;

	return;
}

#undef NRANSI

