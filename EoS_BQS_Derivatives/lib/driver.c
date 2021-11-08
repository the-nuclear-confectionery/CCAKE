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


	/* Vectors for coefficients at low T. */
	// Order 0
	CHI000PAR_ABC=vector(1,3);
	// Order 2
	// Diagonal 
	CHI200PAR_ABC=vector(1,3); CHI020PAR_ABC=vector(1,3); CHI002PAR_ABC=vector(1,3); 
	// Mixed
	CHI110PAR_ABC=vector(1,3);	CHI101PAR_ABC=vector(1,3);	CHI011PAR_ABC=vector(1,3);	
	// Order 4
	// Diagonal
	CHI400PAR_ABC=vector(1,3);	CHI040PAR_ABC=vector(1,3);	CHI004PAR_ABC=vector(1,3);
	// Mixed 31
	CHI310PAR_ABC=vector(1,3);	CHI301PAR_ABC=vector(1,3);	CHI031PAR_ABC=vector(1,3);	
	// Mixed 13
	CHI130PAR_ABC=vector(1,3);	CHI103PAR_ABC=vector(1,3); CHI013PAR_ABC=vector(1,3);	
   // Mixed 22
   CHI220PAR_ABC=vector(1,3); CHI202PAR_ABC=vector(1,3); CHI022PAR_ABC=vector(1,3);	
   // Mied 112
   CHI211PAR_ABC=vector(1,3);	CHI121PAR_ABC=vector(1,3); CHI112PAR_ABC=vector(1,3);


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


	T_min_matching = 70.0;
  	set_lowT_parameters(CHI000PAR, CHI000PAR_ABC);
  	set_lowT_Mod_parameters(CHI200PAR, CHI200PAR_ABC);
  	set_lowT_parameters(CHI020PAR, CHI020PAR_ABC);
  	set_lowT_parameters(CHI002PAR, CHI002PAR_ABC);
  	set_lowT_parameters(CHI110PAR, CHI110PAR_ABC);
  	set_lowT_parameters(CHI101PAR, CHI101PAR_ABC);
  	set_lowT_parameters(CHI011PAR, CHI011PAR_ABC);
  	set_lowT_parameters(CHI400PAR, CHI400PAR_ABC);
  	set_lowT_parameters(CHI040PAR, CHI040PAR_ABC);
  	set_lowT_parameters(CHI004PAR, CHI004PAR_ABC);
  	set_lowT_parameters(CHI310PAR, CHI310PAR_ABC);
  	set_lowT_parameters(CHI301PAR, CHI301PAR_ABC);
  	set_lowT_parameters(CHI031PAR, CHI031PAR_ABC);
  	set_lowT_parameters(CHI130PAR, CHI130PAR_ABC);
  	set_lowT_parameters(CHI103PAR, CHI103PAR_ABC);
  	set_lowT_parameters(CHI013PAR, CHI013PAR_ABC);
  	set_lowT_parameters(CHI220PAR, CHI220PAR_ABC);
  	set_lowT_parameters(CHI202PAR, CHI202PAR_ABC);
  	set_lowT_parameters(CHI022PAR, CHI022PAR_ABC);
  	set_lowT_parameters(CHI211PAR, CHI211PAR_ABC);
  	set_lowT_parameters(CHI121PAR, CHI121PAR_ABC);
  	set_lowT_parameters(CHI112PAR, CHI112PAR_ABC);

  	
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Initialization completed in %lf seconds.\n", cpu_time_used);	

	return;
}




/* The main body of the program. */
void initialize_thermodynamics(const char parameters_filename[])
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


	/* Vectors for coefficients at low T. */
	// Order 0
	CHI000PAR_ABC=vector(1,3);
	// Order 2
	// Diagonal 
	CHI200PAR_ABC=vector(1,3); CHI020PAR_ABC=vector(1,3); CHI002PAR_ABC=vector(1,3); 
	// Mixed
	CHI110PAR_ABC=vector(1,3);	CHI101PAR_ABC=vector(1,3);	CHI011PAR_ABC=vector(1,3);	
	// Order 4
	// Diagonal
	CHI400PAR_ABC=vector(1,3);	CHI040PAR_ABC=vector(1,3);	CHI004PAR_ABC=vector(1,3);
	// Mixed 31
	CHI310PAR_ABC=vector(1,3);	CHI301PAR_ABC=vector(1,3);	CHI031PAR_ABC=vector(1,3);	
	// Mixed 13
	CHI130PAR_ABC=vector(1,3);	CHI103PAR_ABC=vector(1,3); CHI013PAR_ABC=vector(1,3);	
   // Mixed 22
   CHI220PAR_ABC=vector(1,3); CHI202PAR_ABC=vector(1,3); CHI022PAR_ABC=vector(1,3);	
   // Mied 112
   CHI211PAR_ABC=vector(1,3);	CHI121PAR_ABC=vector(1,3); CHI112PAR_ABC=vector(1,3);


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


	T_min_matching = 70.0;
  	set_lowT_parameters(CHI000PAR, CHI000PAR_ABC);
  	set_lowT_Mod_parameters(CHI200PAR, CHI200PAR_ABC);
  	set_lowT_parameters(CHI020PAR, CHI020PAR_ABC);
  	set_lowT_parameters(CHI002PAR, CHI002PAR_ABC);
  	set_lowT_parameters(CHI110PAR, CHI110PAR_ABC);
  	set_lowT_parameters(CHI101PAR, CHI101PAR_ABC);
  	set_lowT_parameters(CHI011PAR, CHI011PAR_ABC);
  	set_lowT_parameters(CHI400PAR, CHI400PAR_ABC);
  	set_lowT_parameters(CHI040PAR, CHI040PAR_ABC);
  	set_lowT_parameters(CHI004PAR, CHI004PAR_ABC);
  	set_lowT_parameters(CHI310PAR, CHI310PAR_ABC);
  	set_lowT_parameters(CHI301PAR, CHI301PAR_ABC);
  	set_lowT_parameters(CHI031PAR, CHI031PAR_ABC);
  	set_lowT_parameters(CHI130PAR, CHI130PAR_ABC);
  	set_lowT_parameters(CHI103PAR, CHI103PAR_ABC);
  	set_lowT_parameters(CHI013PAR, CHI013PAR_ABC);
  	set_lowT_parameters(CHI220PAR, CHI220PAR_ABC);
  	set_lowT_parameters(CHI202PAR, CHI202PAR_ABC);
  	set_lowT_parameters(CHI022PAR, CHI022PAR_ABC);
  	set_lowT_parameters(CHI211PAR, CHI211PAR_ABC);
  	set_lowT_parameters(CHI121PAR, CHI121PAR_ABC);
  	set_lowT_parameters(CHI112PAR, CHI112PAR_ABC);

  	
	end = clock();
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Initialization completed in %lf seconds.\n", cpu_time_used);	

	return;
}


/*
void get_densities(double point[], double densities[])
{
	const double hbarc3 = 197.3*197.3*197.3;
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
*/



//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// S <<-->> Q HAVE BEEN REVERSED!!!!!!!
void get_eBSQ_densities(double point[], double densities[])
{
	const double Tsol = point[0], muBsol = point[1], muSsol = point[2], muQsol = point[3];
	const double Tsol3_by_hc3 = Tsol*Tsol*Tsol/(197.3*197.3*197.3);
	double POut = Tsol*Tsol3_by_hc3*PressTaylor(Tsol, muBsol, muQsol, muSsol);
	double sOut = Tsol3_by_hc3*EntrTaylor(Tsol, muBsol, muQsol, muSsol);
	double BOut = Tsol3_by_hc3*BarDensTaylor(Tsol, muBsol, muQsol, muSsol);
	double SOut = Tsol3_by_hc3*StrDensTaylor(Tsol, muBsol, muQsol, muSsol);
	double QOut = Tsol3_by_hc3*ChDensTaylor(Tsol, muBsol, muQsol, muSsol);
	double eOut = sOut*Tsol - POut + muBsol*BOut + muQsol*QOut + muSsol*SOut;
	densities[0] = eOut;
	densities[1] = BOut;
	densities[2] = SOut;
	densities[3] = QOut;
}

void get_sBSQ_densities(double point[], double densities[])
{
	const double Tsol = point[0], muBsol = point[1], muSsol = point[2], muQsol = point[3];
	const double Tsol3_by_hc3 = Tsol*Tsol*Tsol/(197.3*197.3*197.3);
	densities[0] = Tsol3_by_hc3*EntrTaylor(Tsol, muBsol, muQsol, muSsol);
	densities[1] = Tsol3_by_hc3*BarDensTaylor(Tsol, muBsol, muQsol, muSsol);
	densities[2] = Tsol3_by_hc3*StrDensTaylor(Tsol, muBsol, muQsol, muSsol);
	densities[3] = Tsol3_by_hc3*ChDensTaylor(Tsol, muBsol, muQsol, muSsol);
}



void get_full_thermo(double point[], double thermodynamics[])
{
	const double Tsol = point[0], muBsol = point[1], muSsol = point[2], muQsol = point[3];
	const double Tsol2_by_hc2 = Tsol*Tsol/(197.3*197.3);
	const double Tsol3_by_hc3 = Tsol*Tsol*Tsol/(197.3*197.3*197.3);
	double POut = Tsol*Tsol3_by_hc3*PressTaylor(Tsol, muBsol, muQsol, muSsol);
	double sOut = Tsol3_by_hc3*EntrTaylor(Tsol, muBsol, muQsol, muSsol);
	double BOut = Tsol3_by_hc3*BarDensTaylor(Tsol, muBsol, muQsol, muSsol);
	double SOut = Tsol3_by_hc3*StrDensTaylor(Tsol, muBsol, muQsol, muSsol);
	double QOut = Tsol3_by_hc3*ChDensTaylor(Tsol, muBsol, muQsol, muSsol);
	double eOut = sOut*Tsol - POut + muBsol*BOut + muQsol*QOut + muSsol*SOut;


	//Thermodynamics
	thermodynamics[0]  = POut / 197.3;
	thermodynamics[1]  = sOut;
	thermodynamics[2]  = BOut;
	thermodynamics[3]  = SOut;
	thermodynamics[4]  = QOut;
	thermodynamics[5]  = eOut / 197.3;
	thermodynamics[6]  = SpSound(Tsol, muBsol, muQsol, muSsol);
				
	//Second Order Derivatives (prefactor converts to physical susceptibilities in fm^-2)
	thermodynamics[7]  = Tsol2_by_hc2*P2B2(Tsol, muBsol, muQsol, muSsol);
	thermodynamics[8]  = Tsol2_by_hc2*P2Q2(Tsol, muBsol, muQsol, muSsol);
	thermodynamics[9]  = Tsol2_by_hc2*P2S2(Tsol, muBsol, muQsol, muSsol);
	
	thermodynamics[10] = Tsol2_by_hc2*P2BQ(Tsol, muBsol, muQsol, muSsol);
	thermodynamics[11] = Tsol2_by_hc2*P2BS(Tsol, muBsol, muQsol, muSsol);
	thermodynamics[12] = Tsol2_by_hc2*P2QS(Tsol, muBsol, muQsol, muSsol);
	
	thermodynamics[13] = Tsol2_by_hc2*P2TB(Tsol, muBsol, muQsol, muSsol);
	thermodynamics[14] = Tsol2_by_hc2*P2TQ(Tsol, muBsol, muQsol, muSsol);
	thermodynamics[15] = Tsol2_by_hc2*P2TS(Tsol, muBsol, muQsol, muSsol);
	thermodynamics[16] = Tsol2_by_hc2*P2T2(Tsol, muBsol, muQsol, muSsol);
}
























////////////////////////////////////////////////////////////////////////////////
// STANDARD ARGUMENT ORDERING BELOW THIS LINE

void STANDARD_get_eBSQ_densities(double point[], double densities[])
{
	const double Tsol = point[0]*197.3, muBsol = point[1]*197.3,
                muQsol = point[2]*197.3, muSsol = point[3]*197.3;
	const double Tsol3_by_hc3 = Tsol*Tsol*Tsol/(197.3*197.3*197.3);
	double POut = Tsol*Tsol3_by_hc3*PressTaylor(Tsol, muBsol, muQsol, muSsol);
	double sOut = Tsol3_by_hc3*EntrTaylor(Tsol, muBsol, muQsol, muSsol);
	double BOut = Tsol3_by_hc3*BarDensTaylor(Tsol, muBsol, muQsol, muSsol);
	double SOut = Tsol3_by_hc3*StrDensTaylor(Tsol, muBsol, muQsol, muSsol);
	double QOut = Tsol3_by_hc3*ChDensTaylor(Tsol, muBsol, muQsol, muSsol);
	double eOut = sOut*Tsol - POut + muBsol*BOut + muQsol*QOut + muSsol*SOut;
	densities[0] = eOut / 197.3;
	densities[1] = BOut;
	densities[2] = SOut;
	densities[3] = QOut;
}

void STANDARD_get_sBSQ_densities(double point[], double densities[])
{
	const double Tsol = point[0]*197.3, muBsol = point[1]*197.3,
                muQsol = point[2]*197.3, muSsol = point[3]*197.3;
	const double Tsol3_by_hc3 = Tsol*Tsol*Tsol/(197.3*197.3*197.3);
	densities[0] = Tsol3_by_hc3*EntrTaylor(Tsol, muBsol, muQsol, muSsol);
	densities[1] = Tsol3_by_hc3*BarDensTaylor(Tsol, muBsol, muQsol, muSsol);
	densities[2] = Tsol3_by_hc3*StrDensTaylor(Tsol, muBsol, muQsol, muSsol);
	densities[3] = Tsol3_by_hc3*ChDensTaylor(Tsol, muBsol, muQsol, muSsol);
}



void STANDARD_get_full_thermo(double point[], double thermodynamics[])
{
	const double Tsol = point[0]*197.3, muBsol = point[1]*197.3,
                muQsol = point[2]*197.3, muSsol = point[3]*197.3;
	const double Tsol2_by_hc2 = Tsol*Tsol/(197.3*197.3);
	const double Tsol3_by_hc3 = Tsol*Tsol*Tsol/(197.3*197.3*197.3);
	double POut = Tsol*Tsol3_by_hc3*PressTaylor(Tsol, muBsol, muQsol, muSsol);
	double sOut = Tsol3_by_hc3*EntrTaylor(Tsol, muBsol, muQsol, muSsol);
	double BOut = Tsol3_by_hc3*BarDensTaylor(Tsol, muBsol, muQsol, muSsol);
	double SOut = Tsol3_by_hc3*StrDensTaylor(Tsol, muBsol, muQsol, muSsol);
	double QOut = Tsol3_by_hc3*ChDensTaylor(Tsol, muBsol, muQsol, muSsol);
	double eOut = sOut*Tsol - POut + muBsol*BOut + muQsol*QOut + muSsol*SOut;


	//Thermodynamics
	thermodynamics[0]  = POut / 197.3;
	thermodynamics[1]  = sOut;
	thermodynamics[2]  = BOut;
	thermodynamics[3]  = SOut;
	thermodynamics[4]  = QOut;
	thermodynamics[5]  = eOut / 197.3;
	thermodynamics[6]  = SpSound(Tsol, muBsol, muQsol, muSsol);
				
	//Second Order Derivatives (prefactor converts to physical susceptibilities in fm^-2)
	thermodynamics[7]  = Tsol2_by_hc2*P2B2(Tsol, muBsol, muQsol, muSsol);
	thermodynamics[8]  = Tsol2_by_hc2*P2Q2(Tsol, muBsol, muQsol, muSsol);
	thermodynamics[9]  = Tsol2_by_hc2*P2S2(Tsol, muBsol, muQsol, muSsol);
	
	thermodynamics[10] = Tsol2_by_hc2*P2BQ(Tsol, muBsol, muQsol, muSsol);
	thermodynamics[11] = Tsol2_by_hc2*P2BS(Tsol, muBsol, muQsol, muSsol);
	thermodynamics[12] = Tsol2_by_hc2*P2QS(Tsol, muBsol, muQsol, muSsol);
	
	thermodynamics[13] = Tsol2_by_hc2*P2TB(Tsol, muBsol, muQsol, muSsol);
	thermodynamics[14] = Tsol2_by_hc2*P2TQ(Tsol, muBsol, muQsol, muSsol);
	thermodynamics[15] = Tsol2_by_hc2*P2TS(Tsol, muBsol, muQsol, muSsol);
	thermodynamics[16] = Tsol2_by_hc2*P2T2(Tsol, muBsol, muQsol, muSsol);
}



#undef NRANSI
