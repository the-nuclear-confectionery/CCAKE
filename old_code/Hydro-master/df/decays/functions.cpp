#ifndef FUNCTIONS_CPP_
#define FUNCTIONS_CPP_


/* functions.c            functions for reso.c main program calculations
*   Josef Sollfrank                  Nov. .98             */

// Commented and adjusted by Evan Frodermann, Aug 2005

#include	<string.h>
#include	<stdio.h>
#include	<math.h>
#include        <stdlib.h>
#include <iostream>
#include <fstream>


#include	"functions.h"
#include	"tools.h"
#include	"reso.h"



using namespace std;

//#include        <gsl/gsl_sf_bessel.h>

#define         PTCHANGE    1.0

// double thermalspectra(double pt, double T, double mass, double mu, double degen);

//*******************************************************************************
// The following arrays are set for specific integration routines and point with
// weights for those routines

static	double	gaus2x[] = { 0.577350269189626 };
static	double	gaus4x[] = {	0.8611363115,	0.3399810435	};
static	double	gaus8x[] = {	0.9602898564,	0.7966664774,
				0.3137066458,	0.3626837833	};
static	double	gaus10x[] = {	0.1488743389,	0.4333953941,
				0.6794095682,	0.8650633666,
						0.97390652	};
static	double	gaus12x[] = {	0.9815606342,	0.9041172563,
				0.7699026741,	0.5873179542,
				0.3678314989,	0.1252334085	};
static	double	gaus16x[] = {
		0.989400934991650,	0.944575023073233,
		0.865631202387832,	0.755404408355003,
		0.617876244402644,	0.458016777657227,
		0.281603550779259,	0.095012509837637	};
static	double	gaus20x[] = {
		0.993128599185094,	0.963971927277913,
		0.912234428251325,	0.839116971822218,
		0.746331906460150,	0.636053680726515,
		0.510867001950827,	0.373706088715419,
		0.227785851141645,	0.076526521133497	};
static	double	gaus48x[] = {
		0.998771007252426118601,	0.993530172266350757548,
		0.984124583722826857745,	0.970591592546247250461,
		0.952987703160430860723,	0.931386690706554333114,
		0.905879136715569672822,	0.876572020274247885906,
		0.843588261624393530711,	0.807066204029442627083,
		0.767159032515740339254,	0.724034130923814654674,
		0.677872379632663905212,	0.628867396776513623995,
		0.577224726083972703818,	0.523160974722233033678,
		0.466902904750958404545,	0.408686481990716729916,
		0.348755886292160738160,	0.287362487355455576736,
		0.224763790394689061225,	0.161222356068891718056,
		0.097004699209462698930,	0.032380170962869362033 };

static	double
gala4x[] = {	0.322547689619,		1.745761101158,
		4.536620296921,		9.395070912301	};

static	double
gala8x[] = {	0.170279632305,		0.903701776799,
		2.251086629866,		4.266700170288,
		7.045905402393,		10.758516010181,
		15.740678641278,	22.863131736889	};
static	double
gala12x[] = {	0.115722117358,		0.611757484515,
		1.512610269776,		2.833751337744,
		4.599227639418,		6.844525453115,
		9.621316842457,		13.006054993306,
		17.116855187462,	22.151090379397,
		28.487967250984,	37.099121044467	};

static	double
gala15x[] = {	0.093307812017,         0.492691740302,
		1.215595412071,         2.269949526204,
		3.667622721751,         5.425336627414,
		7.565916226613,        10.120228568019,
	       13.130282482176,        16.654407708330,
               20.776478899449,        25.623894226729,
               31.407519169754,        38.530683306486,
	       48.026085572686	};

static	double
gala48x[] = { 2.9811235829960e-02,   0.15710799061788,
    0.38626503757646,   0.71757469411697,
     1.1513938340264,    1.6881858234190,
     2.3285270066532,    3.0731108616526,
     3.9227524130465,    4.8783933559213,
     5.9411080546246,    7.1121105358907,
     8.3927625990912,    9.7845831846873,
     11.289259168010,    12.908657778286,
     14.644840883210,    16.500081428965,
     18.476882386874,    20.577998634022,
     22.806462290521,    25.165612156439,
     27.659128044481,    30.291071001009,
     33.065930662499,    35.988681327479,
     39.064848764198,    42.300590362903,
     45.702792038511,    49.279186382837,
     53.038498087817,    56.990624814804,
     61.146864786140,    65.520206929019,
     70.125706236113,    74.980977518911,
     80.106857350324,    85.528311116034,
     91.275707993668,    97.386667713582,
    103.908833357176,    110.90422088498,
     118.45642504628,    126.68342576889,
     135.76258957786,    145.98643270946,
     157.91561202298,    172.99632814856 };
//*********************************************************************************




/**************************************************************************
*									  *
*   readin() 								  *
*									  *
*   reads in the particle data file and stores the datas in the arrays	  *
*   particle.* and decay.* and fills up the rest of data (antibaryons)	  *
***********************************************************************   *
**************************************************************************/

void   readin(char * filename, int *particlemax, int *decaymax)
{
	 int i=0, j=0, k, h;
	 FILE	*dat;
	 double dummy1;

	 for(k=0;k<MAXINTV;k++)
	   partid[k] = -1;

	 dat = fopen(filename,"r");
	 if(dat == NULL){
	   printf(" NO file: %s  available ! \n", filename);
	   printf(" GOOD BYE AND HAVE A NICE DAY! \n");
	   exit(0);
	}
	 // Read in the particle data from the specified resonance table
	 // Save the data is the structure particle[pn]


	while(fscanf(dat,"%i%s%lf%lf%i%i%i%i%i%lf%i%i", &particle[i].monval,
	       particle[i].name, &particle[i].mass, &particle[i].width,
	       &particle[i].gspin, &particle[i].baryon,
	       &particle[i].strange, &particle[i].charm,
	       &particle[i].bottom, &particle[i].gisospin,
	       &particle[i].charge, &particle[i].decays)==12)
	  {



	   // skips the readin of anti-baryons from resoweak since they are added with the baryon
		if (particle[i].baryon<0){

			for(k=0;k<particle[i].decays;k++)
	     		{
	       h=fscanf(dat,"%i%i%lf%i%i%i%i%i",
			&decayL[j].reso, &decayL[j].numpart, &decayL[j].branch,
			&decayL[j].part[0], &decayL[j].part[1], &decayL[j].part[2],
			&decayL[j].part[3], &decayL[j].part[4]);
			//if (i==100) cout << particle[i].name << " " << particle[i].decays << " " << j << " " << decay[j].numpart << endl;
			//if (decay[j].reso==-32214) cout << j << " " <<  decay[j].numpart << " " << decay[j].branch << " " << decay[j].part[0] << " " << decay[j].part[1] << " " << decay[j].part[2] << endl;
			j++;


			}
			continue;

		}


	    partid[MHALF + particle[i].monval] = i;
	    particle[i].stable = 0;


	    /* read in the decays */
	    // These decays are saved in a seperate data set, decay[i].

	   for(k=0;k<particle[i].decays;k++)
	     {
	       h=fscanf(dat,"%i%i%lf%i%i%i%i%i",
			&decayL[j].reso, &decayL[j].numpart, &decayL[j].branch,
			&decayL[j].part[0], &decayL[j].part[1], &decayL[j].part[2],
			&decayL[j].part[3], &decayL[j].part[4]);

			//if (i==98) cout << particle[i].name << " " << particle[i].decays << " " << j << " " << decay[j].numpart << endl;

			if (decayL[j].numpart>=5) cout << decayL[j].reso << " " << decayL[j].numpart << endl;


	        if (h != 8) {
	          printf("Error in scanf decay \n");
                  printf(" GOOD BYE AND HAVE A NICE DAY! \n");
	          exit(0);
	        }
		if (decayL[j].numpart == 1) { particle[i].stable = 1;
		cout << "stable: " <<  particle[i].monval << " " << particle[i].name << endl;}
		j++; // Add one to the decay counting variable "j"
	    }
	   //printf("\n");

	   /* setting of additional parameters */

	    if (particle[i].baryon == 1)
	      {
		i++;// If the particle is a baryon, add a particle for the anti-baryon
		    // Add one to the counting variable "i" for the number of particles for the antibaryon
		particle[i].monval = -particle[i-1].monval;
		strcpy(particle[i].name, "  Anti-");
		strncat(particle[i].name, particle[i-1].name, 18);
		particle[i].mass     =  particle[i-1].mass;
		particle[i].width    =  particle[i-1].width;
		particle[i].gspin    =  particle[i-1].gspin;
		particle[i].baryon   = -particle[i-1].baryon;
		particle[i].strange  = -particle[i-1].strange;
		particle[i].charm    = -particle[i-1].charm;
		particle[i].bottom   = -particle[i-1].bottom;
		particle[i].gisospin =  particle[i-1].gisospin;
		particle[i].charge   = -particle[i-1].charge;
		particle[i].decays   = particle[i-1].decays;
	 	partid[MHALF + particle[i].monval] = i;
	        particle[i].stable =  particle[i-1].stable;


	    }
	    i++; // Add one to the counting variable "i" for the meson/baryon
	  }
	fclose(dat);

	//printf("last particle: %i %s \n",particle[i-1].monval,
        //                  particle[i-1].name);
	//printf("# particle %5i; # decays %5i; \n\n",i,j);

	*particlemax = i;   // Set the maxparticle variable to be the value of "i"
	*decaymax  = j;     // Set the maxdecays variable to be the value of "j"

	cout << "decay max " << j << endl;

        if((*particlemax) > NUMPARTICLE){
          printf("Array for particles to small!!\n");
          printf(" GOOD BYE AND HAVE A NICE DAY! \n");
	  exit(0);
	}
      }

//*******************************************************************************************************
// This function reads in the dN/d3p spectra calculated by spectra via the HYDRO output.  Set the pt
// points with the gaussian values given above.  These are the same values given to the spectra in the
// azspectra0p0 program. The default filename is phipspectra.dat

void   readspec( string specfile , char resofile [], int *particlemax, int *decaymax,int nphi,int npt,int ptmax)

/*       char  infile[FILEDIM];*/
/*       char  specfile[FILEDIM];*/
/*       int *particlemax, *decaymax;*/
       {

       // FILE    *dat;
        FILE    *spec;

        int i, j, k, pn, npa = 0;
        int mon;



//	dat = fopen(infile,"r");
//        if(dat == NULL){
//          printf(" NO file: %s  available ! \n", infile);
//          printf(" GOOD BYE AND HAVE A NICE DAY! \n");
//          exit(0);
//	}

//     setup gamma
	 npa++;
        pn = partid[MHALF + 22];
	switch (npt) {
	    case 4:
              for(i=0;i<4;i++)
                particle[pn].pt[i] =  ptmax*gala4x[i]/gala4x[3];
            break;
	    case 8:
              for(i=0;i<8;i++)
                particle[pn].pt[i] =  ptmax*gala8x[i]/ gala8x[7];
            break;
	    case 12:
              for(i=0;i<12;i++)
                particle[pn].pt[i] =  ptmax*gala12x[i]/gala12x[12];
            break;
	    case 15:
              for(i=0;i<15;i++)
                particle[pn].pt[i] =  ptmax*gala15x[i]/gala15x[14];
            break;
	    case 48:
              for(i=0;i<48;i++)
                particle[pn].pt[i] =  ptmax*gala48x[i]/gala48x[47];
            break;
            default:
              printf(" No abscissas for npt = %i !\n",npt);
              printf(" GOOD BYE AND HAVE A NICE DAY! \n");
              exit(0);
	  }

	 particle[pn].nphi = nphi;
          particle[pn].npt = npt;
          for(k=0;k<nphi;k++){
            for(j=0;j<npt;j++) particle[pn].dNdptdphi[j][k]=0;
		}






//      setup gamma




	spec = fopen(specfile.c_str(),"r");
        if(spec == NULL){
          printf(" NO file: %s  available ! \n", specfile.c_str());
          printf(" GOOD BYE AND HAVE A NICE DAY! \n");
          exit(0);
	}
        //double dum=fscanf(dat,"%s",resofile);
        readin(resofile, particlemax, decaymax);

        while( fscanf(spec,"%i",&mon) == 1){
       // cout << mon << endl;


          pn = partid[MHALF + mon];
         // cout << pn << " " << mon << endl;
          //printf(" %i %i %lf %i %i %lf %lf \n",npa, mon,  slope, npt, nphi,
          //        dum1, dum2);
          if(pn == -1){
            printf(" particle %i not in reso table ! \n",mon);
            continue;
	  }

        npa++;

          if(npt > NPT){
            printf(" NPT = %i array to small !\n", npt);
            printf(" GOOD BYE AND HAVE A NICE DAY! \n");
            exit(0);
	  }
          if(nphi > NPHI){
            printf(" NPHI = %i array to small !\n", nphi);
            printf(" GOOD BYE AND HAVE A NICE DAY! \n");
            exit(0);
	  }
          switch (npt) {
	    case 4:
              for(i=0;i<4;i++)
                particle[pn].pt[i] =  ptmax*gala4x[i]/gala4x[3];
            break;
	    case 8:
              for(i=0;i<8;i++)
                particle[pn].pt[i] =  ptmax*gala8x[i]/ gala8x[7];
            break;
	    case 12:
              for(i=0;i<12;i++)
                particle[pn].pt[i] =  ptmax*gala12x[i]/gala12x[12];
            break;
	    case 15:
              for(i=0;i<15;i++)
                particle[pn].pt[i] =  ptmax*gala15x[i]/gala15x[14];
            break;
	    case 48:
              for(i=0;i<48;i++)
                particle[pn].pt[i] =  ptmax*gala48x[i]/gala48x[47];
            break;
            default:
              printf(" No abscissas for npt = %i !\n",npt);
              printf(" GOOD BYE AND HAVE A NICE DAY! \n");
              exit(0);
	  }
          switch (nphi) {
	    case 2:
              for(i=0;i<1;i++){
                PHI[1-i] = PI*(gaus2x[i] + 1.0);
                PHI[i] = PI*(1.0 - gaus2x[i]);
	      }
            break;
	    case 4:
              for(i=0;i<2;i++){
                PHI[3-i] = PI*(gaus4x[i] + 1.0);
                PHI[i] = PI*(1.0 - gaus4x[i]);
	      }
            break;
	    case 8:
              for(i=0;i<4;i++){
                PHI[7-i] = PI*(gaus8x[i] + 1.0);
                PHI[i] = PI*(1.0 - gaus8x[i]);
	      }
            break;
	    case 10:
              for(i=0;i<5;i++){
                PHI[9-i] = PI*(gaus10x[i] + 1.0);
                PHI[i] =PI*(1.0 - gaus10x[i]);
	      }
            break;
	    case 12:
              for(i=0;i<6;i++){
                PHI[11-i] = PI*(gaus12x[i] + 1.0);
                PHI[i] = PI*(1.0 - gaus12x[i]);
	      }
            break;
	    case 16:
              for(i=0;i<8;i++){
                PHI[15-i] =PI*(gaus16x[i] + 1.0);
                PHI[i] = PI*(1.0 - gaus16x[i]);
	      }
            break;
	    case 20:
              for(i=0;i<10;i++){
                PHI[19-i] = PI*(gaus20x[i] + 1.0);
                PHI[i] = PI*(1.0 - gaus20x[i]);
	      }
            break;
	    case 48:
              for(i=0;i<24;i++){
                PHI[47-i] = PI*(gaus48x[i] + 1.0);
                PHI[i] = PI*(1.0 - gaus48x[i]);
	      }
            break;
            default:
              printf(" No abscissas for nPhi = %i !\n",nphi);
              printf(" GOOD BYE AND HAVE A NICE DAY! \n");
              exit(0);
	  }
	  if (particle[pn].baryon==1) {
	       int  pnant = partid[MHALF - mon];
	  	npa++;
	  	switch (npt) {
			    case 4:
			      for(i=0;i<4;i++)
				particle[pnant].pt[i] =  ptmax*gala4x[i]/gala4x[3];
			    break;
			    case 8:
			      for(i=0;i<8;i++)
				particle[pnant].pt[i] =  ptmax*gala8x[i]/ gala8x[7];
			    break;
			    case 12:
			      for(i=0;i<12;i++)
				particle[pnant].pt[i] =  ptmax*gala12x[i]/gala12x[12];
			    break;
			    case 15:
			      for(i=0;i<15;i++)
				particle[pnant].pt[i] =  ptmax*gala15x[i]/gala15x[14];
			    break;
			    case 48:
			      for(i=0;i<48;i++)
				particle[pnant].pt[i] =  ptmax*gala48x[i]/gala48x[47];
			    break;
			    default:
			      printf(" No abscissas for npt = %i !\n",npt);
			      printf(" GOOD BYE AND HAVE A NICE DAY! \n");
			      exit(0);
			  }

	  }

         // particle[pn].slope[npt] = slope;
          particle[pn].nphi = nphi;
          particle[pn].npt = npt;
          for(k=0;k<nphi;k++){
            for(j=0;j<npt;j++){
              double dum=fscanf(spec,"%lf",&particle[pn].dNdptdphi[j][k] );
           //   if ((isnan(particle[pn].dNdptdphi[j][k])==1)||(particle[pn].dNdptdphi[j][k]<0)) particle[pn].dNdptdphi[j][k]=0;
          	if (isnan(particle[pn].dNdptdphi[j][k])==1) particle[pn].dNdptdphi[j][k]=0;
              	if (particle[pn].baryon==1)  {
              	 int  pnant = partid[MHALF - mon];
              	particle[pnant].dNdptdphi[j][k]=particle[pn].dNdptdphi[j][k];
              	particle[pnant].nphi = nphi;
          	particle[pnant].npt = npt;
              	}

	    }
	  }

	//if (particle[pn].monval==213) cout << particle[pn].dNdptdphi[0][0] << endl;

	}
       // fclose(dat);
        fclose(spec);
//        if(npa > 0)
//          printf(" Successful read in of %5i spectra !\n",npa);
      }

//***********************************************************************************************************
// After calculating the spectra in the decay routines, this routine writes that data to file.  The spectra is written
// to spec_###.dat in block format for the pt/phi dependence.  The pt values for each point are saved to specPT_###.dat.
// The ### corresponds to the monte carlo value number assigned the resonance table.  The phi data was already stored in
// angle.dat. (by default)

void   writespec(int particlemax, string totout)

/*     int particlemax;*/
/*     char outdir[FILEDIM];*/

{





  char p[SSL];
  int i, j, k;


  ofstream OUT;
  OUT.open(totout.c_str());
  for(i=0;i<particlemax;i++)// Cycle through the particles
    {

      if(particle[i].stable == 1||particle[i].monval==213||particle[i].monval==-213) //Only print out the stable particles
	{


	  convei(particle[i].monval,p);


	if ((particle[i].npt>0)&&(particle[i].nphi>0)) OUT <<  p << endl;
	 for(k=0;k<particle[i].nphi;k++)//Print out the desired data.
	    {
	 for(j=0;j<particle[i].npt;j++)
	{
	 	OUT <<  particle[i].dNdptdphi[j][k] << " " ;
	}
		OUT << endl;
	}
//
//	 if (particle[i].baryon == 1){
//
//
//	 	convei(particle[i+1].monval,p);
//
//
//		if ((particle[i].npt>0)&&(particle[i].nphi>0)) OUT <<  p << endl;
//	 	for(k=0;k<particle[i].nphi;k++)//Print out the desired data.
//	    	{
//	 	for(j=0;j<particle[i].npt;j++)
//		{
//	 	OUT <<  particle[i].dNdptdphi[j][k] << " " ;
//		}
//		OUT << endl;
//		}
//
//	  }
	}

    }
    OUT.close();
    printf(" Printed Spectra \n");

//    /// print off pt and phi calc points
//    ofstream OUT2;
//    strcpy(filename2,outdir);
//    strcat(filename2, "/ptphipoints.dat");
//  OUT2.open(filename2);
//  OUT2 << "pt: " <<  particle[1].npt << endl;
//  for(j=0;j<particle[1].npt;j++) OUT2 << particle[1].pt[j] << " " ;
//  OUT2 << endl ;
//  OUT2 << "phi: " << particle[1].nphi << endl;
//  for(k=0;k<particle[1].nphi;k++) OUT2 << PHI[k] << " " ;
//  OUT2 << endl ;
//
//
//    OUT2.close();
//    /// print off pt and phi calc points






}





/*************************************************
*
*	Edndp3
*
*
**************************************************/
// This function interpolates the needed spectra for a given pt and phi.

double	Edndp3(double yr, double  ptr, double phirin,int  res_num)

//				/* supersedes during test the right one */
//	double	yr;		/* y  of resonance */
//	double	ptr;		/* pt of resonance */
//	double	phirin;		/* phi angle  of resonance */
//	int	res_num;	/* Montecarlo number of resonance 	*/

	{
	  double	phir, val;
	  double        f1, f2;
	  int     	pn, npt, nphi;



          if(phirin < 0.0){
             printf("ERROR: phir %15.8le < 0 !!! \n", phirin);exit(0);}
          if(phirin > 2.0*PI){
               printf("ERROR: phir %15.8le > 2PI !!! \n", phirin);exit(0);}
         phir= phirin;

  	  pn = partid[MHALF + res_num];

          nphi = 1;
          while((phir > PHI[nphi])&&(nphi<(particle[pn].nphi-1))) nphi++;
          npt = 1;
          while((ptr > particle[pn].pt[npt]) && (npt<(particle[pn].npt - 1))) npt++;


        /* phi interpolation */

         f1 = lin_int(PHI[nphi-1], PHI[nphi],
                      particle[pn].dNdptdphi[npt-1][nphi-1],
                      particle[pn].dNdptdphi[npt-1][nphi], phir);
         f2 = lin_int(PHI[nphi-1], PHI[nphi],
                      particle[pn].dNdptdphi[npt][nphi-1],
                      particle[pn].dNdptdphi[npt][nphi], phir);

         if (f1<0.) f1=0.; // security: if for some reason we got a negative number of particles (happened in the viscous code at large eta sometimes)
  if (f2<0.) f2=0.;

  	if (isinf(f1)) cout << "f1 infinity" << endl;
  	if (isinf(f2)) cout << "f2 infinity" << endl;


         if(ptr > PTCHANGE && f1!=0 && f2!=0){
           f1 = log(f1);
           f2 = log(f2);
	 }



         val = lin_int(particle[pn].pt[npt-1],particle[pn].pt[npt],
                       f1, f2, ptr);
         if(ptr > PTCHANGE && f1!=0 && f2!=0)
           val = exp(val);

         if (ptr>particle[pn].pt[particle[pn].npt-1]) val=0;

         return val;
	}



//*****************************************************************************************************
/* // Requires either a) A definition of the K bessel function
   // or b) the GSL (Gnu Scientific Library) installed.  Add the flags -lgsl -lgslcblas to the compilation
   // to use the GSL library.

// Thermal function that simulates most spectra given a finte temperature.  Replaced
// by the output of hydro in the real simulation

double thermalspectra(double pt, double T, double mass, double mu, double degen)
{

  double spect = 0.0;
  double mT = sqrt(pt*pt + mass * mass);
  spect = degen/(4*M_PI*M_PI)*exp(mu/T)*mT *gsl_sf_bessel_K1(mT/T);

  return spect;
}
*/

#endif
