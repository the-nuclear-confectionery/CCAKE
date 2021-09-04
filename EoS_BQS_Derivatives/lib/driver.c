/* 
	Copyright (c) 2019, Paolo Parotto and Jamie Stafford, 
	Department of Physics, University of Houston, Houston, TX 77204, US.
*/

/*
	This file produces a Taylor expanded EoS using Lattice QCD coefficients at muB=0.
*/ 

#define NRANSI

#include "driver.h"

/* Import standard libraries. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include "../include/nrD.h"
#include "../include/nrutilD.h"
#include "../include/export_to_hdf.h"
#include "../include/nm_solver.h"

/* Import additional library for many variables used in the program. */
#include "../include/Variables.h"
#include "../include/Functions.h"

/* Strangeness neutrality. */
#define NSN 2

/* Time variables for timing. */
clock_t start, end;
double cpu_time_used;

/* Functions whose zeroes we are seeking for strangeness neutrality. */
void funcvSN(int n,double x[],double f[]){
	f[1] = StrDensTaylor(Tval,muBval,x[1],x[2]);
	f[2] = ChDensTaylor(Tval,muBval,x[1],x[2]) - 0.4*BarDensTaylor(Tval,muBval,x[1],x[2]);
} 

double mapf(double x, double a, double b)
{
	return log(fabs((x-a)/(b-x)));
}

double imapf(double t, double a, double b)
{
	return a + (b-a)/(1.0+exp(-t));
}

void get_seed_points(double a, double b, int n, double result[])
{
	for (int ii = 1; ii < n-1; ii+=1) result[ii] = mapf(a + (b-a)*ii/(n-1.0), a, b);
	result[0] = mapf(a+1e-10, a, b);
	result[n-1] = mapf(b-1e-10, a, b);
}

/* The main body of the program. */
void initialize(const char parameters_filename[])
{
	char buff[FILENAME_MAX];

	/* Define time variables to measure time elapsed creating the tables.*/
	time_t start, stop;

	/* All vectors and matrices are initialized. */
	/* Vectors for coefficients. */
	// Order 0
	CHI000PAR=vector(1,21);
	// Order 2
	// Diagonal 
	CHI200PAR=vector(1,21); CHI020PAR=vector(1,21); CHI002PAR=vector(1,21); 
	// Mixed
	CHI110PAR=vector(1,21);	CHI101PAR=vector(1,21);	CHI011PAR=vector(1,21);	
	// Order 4
	// Diagonal
	CHI400PAR=vector(1,21);	CHI040PAR=vector(1,21);	CHI004PAR=vector(1,21);
	// Mixed 31
	CHI310PAR=vector(1,21);	CHI301PAR=vector(1,21);	CHI031PAR=vector(1,21);	
	// Mixed 13
	CHI130PAR=vector(1,21);	CHI103PAR=vector(1,21); CHI013PAR=vector(1,21);	
   // Mixed 22
   CHI220PAR=vector(1,21); CHI202PAR=vector(1,21); CHI022PAR=vector(1,21);	
   // Mied 112
   CHI211PAR=vector(1,21);	CHI121PAR=vector(1,21); CHI112PAR=vector(1,21);

	/* Matrix for coefficients: 22 coefficients, 21 parameters each. */
	parMatrix=matrix(1,22,1,21);
	
	/* For strangeness neutrality */
	xSN=vector(1,NSN);
 	fSN=vector(1,NSN);
	
	/* Assign the name of the main folder where the program lives and the files we wish to import are located.*/
	getcwd(buff,FILENAME_MAX);
	printf("Current working directory is: \n\n%s \n\n",buff);

   /* Start clock for importing lists */
  	start = clock();

	/* Parametrization parameters are read from the user-input file, and saved. */
	FILE *ParametersIn = fopen(parameters_filename, "r");
  	if (ParametersIn == 0){
  		fprintf(stderr,"failed to open paremeters file\n");
  		exit(1);
  	}
  	for(i=1;fscanf(ParametersIn,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &A0, &A1, &A2, &A3, &A4, &A5, &A6, &A7, 
  	                                                                                &A8, &A9, &B0, &B1, &B2, &B3, &B4, &B5, &B6, &B7, &B8, &B9, &C0) !=EOF; i++){
       parMatrix[i][1] = A0;  	    parMatrix[i][2] = A1;  	    parMatrix[i][3] = A2;  	    parMatrix[i][4] = A3;
  	    parMatrix[i][5] = A4;  	    parMatrix[i][6] = A5;  	    parMatrix[i][7] = A6;  	    parMatrix[i][8] = A7;
  	    parMatrix[i][9] = A8;  	    parMatrix[i][10] = A9;
  	    parMatrix[i][11] = B0;  	    parMatrix[i][12] = B1;  	    parMatrix[i][13] = B2;  	    parMatrix[i][14] = B3;
  	    parMatrix[i][15] = B4;  	    parMatrix[i][16] = B5;  	    parMatrix[i][17] = B6;  	    parMatrix[i][18] = B7;
  	    parMatrix[i][19] = B8;  	    parMatrix[i][20] = B9;
  	    parMatrix[i][21] = C0;
  	    
  	    //printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",A0,A1,A2,A3,A4,A5,A6,A7,B0,B1,B2,B3,B4,B5,B6,B7,C0);
  	}
  	fclose(ParametersIn);	
  	
 
	// Vectors for coefficients.   
	// All vectors are filled with parameters. 
	// Order 0
  	for(i=1;i<=21;i++) CHI000PAR[i] = parMatrix[1][i];
  	// Order 2
	// Diagonal 
  	for(i=1;i<=21;i++) CHI200PAR[i] = parMatrix[2][i];
  	for(i=1;i<=21;i++) CHI020PAR[i] = parMatrix[3][i];
	for(i=1;i<=21;i++) CHI002PAR[i] = parMatrix[4][i];
	// Mixed
	for(i=1;i<=21;i++) CHI110PAR[i] = parMatrix[5][i];
  	for(i=1;i<=21;i++) CHI101PAR[i] = parMatrix[6][i];
  	for(i=1;i<=21;i++) CHI011PAR[i] = parMatrix[7][i];
	// Order 4
   	// Diagonal	
   	for(i=1;i<=21;i++) CHI400PAR[i] = parMatrix[8][i];
  	for(i=1;i<=21;i++) CHI040PAR[i] = parMatrix[9][i];
  	for(i=1;i<=21;i++) CHI004PAR[i] = parMatrix[10][i];
  	// Mixed 31
  	for(i=1;i<=21;i++) CHI310PAR[i] = parMatrix[11][i];
  	for(i=1;i<=21;i++) CHI301PAR[i] = parMatrix[12][i];
  	for(i=1;i<=21;i++) CHI031PAR[i] = parMatrix[13][i];
  	// Mixed 13
  	for(i=1;i<=21;i++) CHI130PAR[i] = parMatrix[14][i];
  	for(i=1;i<=21;i++) CHI103PAR[i] = parMatrix[15][i];
  	for(i=1;i<=21;i++) CHI013PAR[i] = parMatrix[16][i];
  	// Mixed 22
  	for(i=1;i<=21;i++) CHI220PAR[i] = parMatrix[17][i];
  	for(i=1;i<=21;i++) CHI202PAR[i] = parMatrix[18][i];
  	for(i=1;i<=21;i++) CHI022PAR[i] = parMatrix[19][i];
  	// Mixed 211
  	for(i=1;i<=21;i++) CHI211PAR[i] = parMatrix[20][i];
  	for(i=1;i<=21;i++) CHI121PAR[i] = parMatrix[21][i];
  	for(i=1;i<=21;i++) CHI112PAR[i] = parMatrix[22][i];
  	

	/*  	
	int run_density_solver = 1;
	if ( run_density_solver )
	for (int irun = 0; irun <= 2; irun+=1)
	{
		double eIn = 973.563, BIn = -0.316059, SIn = 0.323859, QIn = 1.06384;	// (MeV,1,1,1)/fm^3
		double densities[4] = {eIn, BIn, SIn, QIn};
		double sols[4] = {158.0, -437.0, -112.0, 437.0};		// MeV
		double minima[4] = {155.0, -450.0, -125.0, 425.0};		// MeV
		double maxima[4] = {160.0, -425.0, -100.0, 450.0};		// MeV

		// find the solution
		if (irun==1) solve(densities, sols);
		else if (irun==2)
		{
			start = clock();

			const int n_seeds_per_dimension = 5;	// includes endpoints
			const double dT_seed = (maxima[0] - minima[0])/(n_seeds_per_dimension-1);
			const double dmuB_seed = (maxima[1] - minima[1])/(n_seeds_per_dimension-1);
			const double dmuS_seed = (maxima[2] - minima[2])/(n_seeds_per_dimension-1);
			const double dmuQ_seed = (maxima[3] - minima[3])/(n_seeds_per_dimension-1);

			double seedT[4], seedmuB[4], seedmuS[4], seedmuQ[4];

			get_seed_points(minima[0], maxima[0], n_seeds_per_dimension, seedT);
			get_seed_points(minima[1], maxima[1], n_seeds_per_dimension, seedmuB);
			get_seed_points(minima[2], maxima[2], n_seeds_per_dimension, seedmuS);
			get_seed_points(minima[3], maxima[3], n_seeds_per_dimension, seedmuQ);
			
			int attempts = 0;
			double Tseed = 0.0, muBseed = 0.0, muSseed = 0.0, muQseed = 0.0;
			for (int ii = 1; ii < n_seeds_per_dimension-1; ii+=1)
			for (int jj = 1; jj < n_seeds_per_dimension-1; jj+=1)
			for (int kk = 1; kk < n_seeds_per_dimension-1; kk+=1)
			for (int ll = 1; ll < n_seeds_per_dimension-1; ll+=1)
			{
				Tseed = seedT[ii]; muBseed = seedmuB[jj]; muSseed = seedmuS[kk]; muQseed = seedmuQ[ll];
				double seeds[4] = {Tseed, muBseed, muSseed, muQseed};
				solve2(densities, sols, minima, maxima, seeds);
				attempts += 1;
				if (sols[0] > 0.0) // success
					goto success;
			}
			success:
				end = clock();
				cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
				printf("Finished calculating in %lf seconds.\n", cpu_time_used);	

				printf("Found solution in %d attempts!\n", attempts);
				printf("Seed point: %lf %lf %lf %lf\n",
						imapf(Tseed, minima[0], maxima[0]),
						imapf(muBseed, minima[1], maxima[1]),
						imapf(muSseed, minima[2], maxima[2]),
						imapf(muQseed, minima[3], maxima[3]));
		}
		double Tsol = sols[0], muBsol = sols[1], muSsol = sols[2], muQsol = sols[3];
		Tval = Tsol; muBval = muBsol; muSval = muSsol; muQval = muQsol;
		i = Tsol; j = muBsol; l = muSsol; k = muQsol;	// Q and S reversed
		
		
		printf("Input:\n");
		printf("eIn = %15.8f\n", eIn);
		printf("BIn = %15.8f\n", BIn);
		printf("SIn = %15.8f\n", SIn);
		printf("QIn = %15.8f\n", QIn);
		printf("Solution:\n");
		printf("Tsol = %15.8f\n", Tsol);
		printf("muBsol = %15.8f\n", muBsol);
		printf("muSsol = %15.8f\n", muSsol);
		printf("muQsol = %15.8f\n", muQsol);
		fflush(stdout);
		
		double POut = Tsol*Tsol*Tsol*Tsol*PressTaylor(Tsol, muBsol, muQsol, muSsol);
		double sOut = Tsol*Tsol*Tsol*EntrTaylor(Tsol, muBsol, muQsol, muSsol);
		double BOut = Tsol*Tsol*Tsol*BarDensTaylor(Tsol, muBsol, muQsol, muSsol);
		double SOut = Tsol*Tsol*Tsol*StrDensTaylor(Tsol, muBsol, muQsol, muSsol);
		double QOut = Tsol*Tsol*Tsol*ChDensTaylor(Tsol, muBsol, muQsol, muSsol);
		POut /= 197.327*197.327*197.327;
		sOut /= 197.327*197.327*197.327;
		BOut /= 197.327*197.327*197.327;
		SOut /= 197.327*197.327*197.327;
		QOut /= 197.327*197.327*197.327;
		double eOut = sOut*Tsol - POut + muBsol*BOut + muQsol*QOut + muSsol*SOut;
		
		printf("Check:\n");
		printf("POut = %15.8f\n", POut);
		printf("sOut = %15.8f\n", sOut);
		printf("eOut = %15.8f\n", eOut);
		printf("BOut = %15.8f\n", BOut);
		printf("SOut = %15.8f\n", SOut);
		printf("QOut = %15.8f\n", QOut);
		fflush(stdout);
		
		//if (irun==1) exit(-1);
		
		
		//Thermodynamics
		PressVal = PressTaylor(i,j,k,l);
		EntrVal = EntrTaylor(i,j,k,l);
		BarDensVal = BarDensTaylor(i,j,k,l);
		StrDensVal = StrDensTaylor(i,j,k,l);
		ChDensVal = ChDensTaylor(i,j,k,l);
		EnerDensVal = EntrVal - PressVal 
				+ muBval/Tval*BarDensVal 
				+ muQval/Tval*ChDensVal 
				+ muSval/Tval*StrDensVal;
		SpSoundVal = SpSound(Tval,muBval,muQval,muSval);
		            
		//Second Order Derivatives
		D2PB2 = P2B2(i,j,k,l);
		D2PQ2 = P2Q2(i,j,k,l);
		D2PS2 = P2S2(i,j,k,l);
		
		D2PBQ = P2BQ(i,j,k,l);
		D2PBS = P2BS(i,j,k,l);
		D2PQS = P2QS(i,j,k,l);
		
		D2PTB = P2TB(i,j,k,l);
		D2PTQ = P2TQ(i,j,k,l);
		D2PTS = P2TS(i,j,k,l);
		D2PT2 = P2T2(i,j,k,l);
		
		printf("vals: %lf  %lf  %lf  %lf  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f\n", Tval, muBval, muQval, muSval, PressVal, EntrVal, 
		        BarDensVal, StrDensVal, ChDensVal, EnerDensVal, SpSoundVal);
		printf("derivs: %lf  %lf  %lf  %lf  %3.12f  %3.12f %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f%  3.12f\n", Tval, muBval, muQval, muSval, D2PB2, D2PQ2, D2PS2, D2PBQ, D2PBS, D2PQS,
				D2PTB, D2PTQ, D2PTS, D2PT2);
		
		printf("********************************************************************************\n\n");
		
		fflush(stdout);
		
//		}
	}
	if ( run_density_solver ) exit(-1);

	chdir(buff);

	*/

	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Initialization completed in %lf seconds.\n", cpu_time_used);	

	
	
	return;
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// S <<-->> Q HAVE BEEN REVERSED!!!!!!!
void get_densities(double point[], double densities[])
{
	const double hbarc3 = 197.327*197.327*197.327;
	const double Tsol = point[0], muBsol = point[1], muSsol = point[2], muQsol = point[3];
	double POut = Tsol*Tsol*Tsol*Tsol*PressTaylor(Tsol, muBsol, muQsol, muSsol)/hbarc3;
	double sOut = Tsol*Tsol*Tsol*EntrTaylor(Tsol, muBsol, muQsol, muSsol)/hbarc3;
	double BOut = Tsol*Tsol*Tsol*BarDensTaylor(Tsol, muBsol, muQsol, muSsol)/hbarc3;
	double SOut = Tsol*Tsol*Tsol*StrDensTaylor(Tsol, muBsol, muQsol, muSsol)/hbarc3;
	double QOut = Tsol*Tsol*Tsol*ChDensTaylor(Tsol, muBsol, muQsol, muSsol)/hbarc3;
	double eOut = sOut*Tsol - POut + muBsol*BOut + muQsol*QOut + muSsol*SOut;
	densities[0] = eOut;
	densities[1] = BOut;
	densities[2] = SOut;
	densities[3] = QOut;
}


#undef NRANSI
