#ifndef DECAY_H_ 
#define DECAY_H_

/* decay.h        header file for decay.c 
*  Sollfrank Josef           Nov .98                      */

void	cal_reso_decays(int maxpart,int  maxdecay,int  bound); 
double dnpir1N (double costh, void *para1);
double dn2ptN (double w2, void *para1);
double dn3ptN (double x, void *para1);
double
Edndp3_2bodyN (double y, double pt,double  phi,double  m1,double  m2, double mr, int res_num);
double
Edndp3_3bodyN (double y, double pt, double phi, double m1, double m2,
	       double m3, double mr, double norm3, int res_num);
void add_reso (int pn, int pnR, int k,int  j);
#endif
