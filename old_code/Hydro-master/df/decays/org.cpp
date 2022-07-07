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
          while((ptr > particle[pn].pt[npt]) &&
                (npt<(particle[pn].npt - 1))) npt++; 
        
        /* phi interpolation */

         f1 = lin_int(PHI[nphi-1], PHI[nphi], 
                      particle[pn].dNdptdphi[npt-1][nphi-1], 
                      particle[pn].dNdptdphi[npt-1][nphi], phir);
         f2 = lin_int(PHI[nphi-1], PHI[nphi], 
                      particle[pn].dNdptdphi[npt][nphi-1], 
                      particle[pn].dNdptdphi[npt][nphi], phir);
                      
         if (f1<0.) f1=0.; // security: if for some reason we got a negative number of particles (happened in the viscous code at large eta sometimes)
  if (f2<0.) f2=0.; 
                      

         if(ptr > PTCHANGE){
           f1 = log(f1); 
           f2 = log(f2);
	 }
         val = lin_int(particle[pn].pt[npt-1],particle[pn].pt[npt], 
                       f1, f2, ptr);
         if(ptr > PTCHANGE)
           val = exp(val);
           
          if (isnan(val1))
         { cout << particle[pn].pt[npt] << " " << PHI[nphi] << " " << res_num << endl;
         
         
         }

/*
         printf(" nphi  %i npt %i \n", nphi,npt);
         printf(" f1  %15.8le %15.8le  \n", f1, f2);
         printf(" phi  %15.8lf %15.8lf  \n", PHI[nphi-1], PHI[nphi]); 
         printf(" pt   %15.8lf %15.8lf  \n", particle[pn].pt[npt-1],particle[pn].pt[npt]);
         printf(" phi  %15.8lf pt %15.8lf    val %15.8lf \n", phir, ptr,val); 
         printf(" phi %15.8le %15.8le \n",particle[pn].dNdptdphi[npt][nphi-1],
                                     particle[pn].dNdptdphi[npt][nphi]);
         printf(" pt  %15.8le %15.8le \n",particle[pn].dNdptdphi[npt-1][nphi-1],
                                     particle[pn].dNdptdphi[npt-1][nphi]);
	   	   	
         exit(0);
*/
         return val;
	}

