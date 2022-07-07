/* functions.h        header file for function.c 
*  Sollfrank Josef           Nov .98                      */
#ifndef FUNCTIONS_H_ 
#define FUNCTIONS_H_

#include <string>
#include <string.h>   // Required by strcpy()
#include <stdlib.h>   // Required by malloc()
#include <iostream>


using namespace std;


void   readin(char * filename, int *particlemax, int *decaymax);
void   readspec( string specfile, char resofile [],int *particlemax, int *decaymax,int nphi, int npt,int ptmax);
void   writespec(int particlemax, string outfile);
double	Edndp3(double yr, double  ptr, double phirin,int  res_num);

#endif
