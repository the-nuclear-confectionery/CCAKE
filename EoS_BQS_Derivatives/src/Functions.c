#include <math.h>
#include <stdio.h>
#define NRANSI
#include "../include/Variables.h"


// ------------- THE FUNCTIONAL FORMS --------------- //
// For all CHIS except CHI200
//double coeff(double *par,double x)
//{
//     return par[21] + (par[1] + par[2]*(154.0/x) + par[3]*(154.0/x)*(154.0/x) + par[4]*(154.0/x)*(154.0/x)*(154.0/x) + par[5]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
//                + par[6]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) + par[7]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
//                + par[8]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) + par[9]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
//                + par[10]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x))
//                  *1.0/(par[11] + par[12]*(154.0/x) + par[13]*(154.0/x)*(154.0/x) + par[14]*(154.0/x)*(154.0/x)*(154.0/x)
//                + par[15]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) + par[16]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
//                + par[17]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) + par[18]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)
//                + par[19]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
//                + par[20]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x));
//}
//
//double coeffprime(double *par,double x)
//{
//    return (- 1.0/x*par[2]*(154.0/x) - 2.0/x*par[3]*(154.0/x)*(154.0/x) - 3.0/x*par[4]*(154.0/x)*(154.0/x)*(154.0/x) - 4.0/x*par[5]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
//                - 5.0/x*par[6]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) - 6.0/x*par[7]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
//                - 7.0/x*par[8]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) - 8.0/x*par[9]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
//                - 9.0/x*par[10]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x))
//                     *1.0/(par[11] + par[12]*(154.0/x) + par[13]*(154.0/x)*(154.0/x) + par[14]*(154.0/x)*(154.0/x)*(154.0/x)
//                + par[15]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) + par[16]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
//                + par[17]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) + par[18]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)
//                + par[19]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
//                + par[20]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x))
//           - (par[1] + par[2]*(154.0/x) + par[3]*(154.0/x)*(154.0/x) + par[4]*(154.0/x)*(154.0/x)*(154.0/x) + par[5]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
//                + par[6]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) + par[7]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
//                + par[8]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) + par[9]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
//                + par[10]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x))
//                    *(- 1.0/x*par[12]*(154.0/x) - 2.0/x*par[13]*(154.0/x)*(154.0/x) - 3.0/x*par[14]*(154.0/x)*(154.0/x)*(154.0/x) - 4.0/x*par[15]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
//                - 5.0/x*par[16]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) - 6.0/x*par[17]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
//                - 7.0/x*par[18]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) - 8.0/x*par[19]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)
//                - 9.0/x*par[20]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x))
//                     *1.0/pow((par[11] + par[12]*(154.0/x) + par[13]*(154.0/x)*(154.0/x) + par[14]*(154.0/x)*(154.0/x)*(154.0/x)
//                + par[15]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) + par[16]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
//                + par[17]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) + par[18]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)
//                + par[19]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x) 
//                + par[20]*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)*(154.0/x)),2);
//}


double coeff(double *par,double x)
{
  const double r = 154.0/x;
  const double r2 = r*r;
  const double r3 = r*r2;
  const double r4 = r*r3;
  const double r5 = r*r4;
  const double r6 = r*r5;
  const double r7 = r*r6;
  const double r8 = r*r7;
  const double r9 = r*r8;
  const double ix = 1.0/x;
     return par[21] + (par[1] + par[2]*r + par[3]*r2 + par[4]*r3 + par[5]*r4 
                + par[6]*r5 + par[7]*r6 + par[8]*r7 + par[9]*r8 + par[10]*r9)
                  *1.0/(par[11] + par[12]*r + par[13]*r2 + par[14]*r3
                + par[15]*r4 + par[16]*r5 + par[17]*r6 + par[18]*r7 + par[19]*r8 + par[20]*r9);
}

double coeffprime(double *par,double x)
{
  const double r = 154.0/x;
  const double r2 = r*r;
  const double r3 = r*r2;
  const double r4 = r*r3;
  const double r5 = r*r4;
  const double r6 = r*r5;
  const double r7 = r*r6;
  const double r8 = r*r7;
  const double r9 = r*r8;
  const double ix = 1.0/x;
    return (- 1.0*ix*par[2]*r - 2.0*ix*par[3]*r2 - 3.0*ix*par[4]*r3 - 4.0*ix*par[5]*r4 
                - 5.0*ix*par[6]*r5 - 6.0*ix*par[7]*r6 
                - 7.0*ix*par[8]*r7 - 8.0*ix*par[9]*r8 
                - 9.0*ix*par[10]*r9)
                     *1.0/(par[11] + par[12]*r + par[13]*r2 + par[14]*r3
                + par[15]*r4 + par[16]*r5 
                + par[17]*r6 + par[18]*r7
                + par[19]*r8 
                + par[20]*r9)
           - (par[1] + par[2]*r + par[3]*r2 + par[4]*r3 + par[5]*r4 
                + par[6]*r5 + par[7]*r6 
                + par[8]*r7 + par[9]*r8 
                + par[10]*r9)
                    *(- 1.0*ix*par[12]*r - 2.0*ix*par[13]*r2 - 3.0*ix*par[14]*r3 - 4.0*ix*par[15]*r4 
                - 5.0*ix*par[16]*r5 - 6.0*ix*par[17]*r6 
                - 7.0*ix*par[18]*r7 - 8.0*ix*par[19]*r8
                - 9.0*ix*par[20]*r9)
                     *1.0/pow((par[11] + par[12]*r + par[13]*r2 + par[14]*r3
                + par[15]*r4 + par[16]*r5 
                + par[17]*r6 + par[18]*r7
                + par[19]*r8 
                + par[20]*r9),2);
}


//double coeffsecond(double *par, double x){
//   return (-2*(-((154.0*par[2])/pow(x,2)) - (2*pow(154.0,2)*par[3])/pow(x,3) - (3*pow(154.0,3)*par[4])/pow(x,4) - (4*pow(154.0,4)*par[5])/pow(x,5) - (5*pow(154.0,5)*par[6])/pow(x,6) - 
//        (6*pow(154.0,6)*par[7])/pow(x,7) - (7*pow(154.0,7)*par[8])/pow(x,8) - (8*pow(154.0,8)*par[9])/pow(x,9) - (9*pow(154.0,9)*par[10])/pow(x,10))*
//      (-((154.0*par[12])/pow(x,2)) - (2*pow(154.0,2)*par[13])/pow(x,3) - (3*pow(154.0,3)*par[14])/pow(x,4) - (4*pow(154.0,4)*par[15])/pow(x,5) - (5*pow(154.0,5)*par[16])/pow(x,6) - 
//        (6*pow(154.0,6)*par[17])/pow(x,7) - (7*pow(154.0,7)*par[18])/pow(x,8) - (8*pow(154.0,8)*par[19])/pow(x,9) - (9*pow(154.0,9)*par[20])/pow(x,10)))/
//    pow(par[11] + (154.0*par[12])/x + (pow(154.0,2)*par[13])/pow(x,2) + (pow(154.0,3)*par[14])/pow(x,3) + (pow(154.0,4)*par[15])/pow(x,4) + (pow(154.0,5)*par[16])/pow(x,5) + 
//      (pow(154.0,6)*par[17])/pow(x,6) + (pow(154.0,7)*par[18])/pow(x,7) + (pow(154.0,8)*par[19])/pow(x,8) + (pow(154.0,9)*par[20])/pow(x,9),2) + 
//   ((2*154.0*par[2])/pow(x,3) + (6*pow(154.0,2)*par[3])/pow(x,4) + (12*pow(154.0,3)*par[4])/pow(x,5) + 
//      (20*pow(154.0,4)*par[5])/pow(x,6) + (30*pow(154.0,5)*par[6])/pow(x,7) + (42*pow(154.0,6)*par[7])/pow(x,8) + (56*pow(154.0,7)*par[8])/pow(x,9) + (72*pow(154.0,8)*par[9])/pow(x,10) + 
//      (90*pow(154.0,9)*par[10])/pow(x,11))/
//    (par[11] + (154.0*par[12])/x + (pow(154.0,2)*par[13])/pow(x,2) + (pow(154.0,3)*par[14])/pow(x,3) + (pow(154.0,4)*par[15])/pow(x,4) + (pow(154.0,5)*par[16])/pow(x,5) + 
//      (pow(154.0,6)*par[17])/pow(x,6) + (pow(154.0,7)*par[18])/pow(x,7) + (pow(154.0,8)*par[19])/pow(x,8) + (pow(154.0,9)*par[20])/pow(x,9)) + 
//   (par[1] + (154.0*par[2])/x + (pow(154.0,2)*par[3])/pow(x,2) + (pow(154.0,3)*par[4])/pow(x,3) + (pow(154.0,4)*par[5])/pow(x,4) + (pow(154.0,5)*par[6])/pow(x,5) + 
//      (pow(154.0,6)*par[7])/pow(x,6) + (pow(154.0,7)*par[8])/pow(x,7) + (pow(154.0,8)*par[9])/pow(x,8) + (pow(154.0,9)*par[10])/pow(x,9))*
//    ((2*pow(-((154.0*par[12])/pow(x,2)) - (2*pow(154.0,2)*par[13])/pow(x,3) - (3*pow(154.0,3)*par[14])/pow(x,4) - (4*pow(154.0,4)*par[15])/pow(x,5) - (5*pow(154.0,5)*par[16])/pow(x,6) - 
//           (6*pow(154.0,6)*par[17])/pow(x,7) - (7*pow(154.0,7)*par[18])/pow(x,8) - (8*pow(154.0,8)*par[19])/pow(x,9) - (9*pow(154.0,9)*par[20])/pow(x,10),2))/
//       pow(par[11] + (154.0*par[12])/x + (pow(154.0,2)*par[13])/pow(x,2) + (pow(154.0,3)*par[14])/pow(x,3) + (pow(154.0,4)*par[15])/pow(x,4) + (pow(154.0,5)*par[16])/pow(x,5) + 
//         (pow(154.0,6)*par[17])/pow(x,6) + (pow(154.0,7)*par[18])/pow(x,7) + (pow(154.0,8)*par[19])/pow(x,8) + (pow(154.0,9)*par[20])/pow(x,9),3) - 
//      ((2*154.0*par[12])/pow(x,3) + (6*pow(154.0,2)*par[13])/pow(x,4) + (12*pow(154.0,3)*par[14])/pow(x,5) + (20*pow(154.0,4)*par[15])/pow(x,6) + (30*pow(154.0,5)*par[16])/pow(x,7) + 
//         (42*pow(154.0,6)*par[17])/pow(x,8) + (56*pow(154.0,7)*par[18])/pow(x,9) + (72*pow(154.0,8)*par[19])/pow(x,10) + (90*pow(154.0,9)*par[20])/pow(x,11))/
//       pow(par[11] + (154.0*par[12])/x + (pow(154.0,2)*par[13])/pow(x,2) + (pow(154.0,3)*par[14])/pow(x,3) + (pow(154.0,4)*par[15])/pow(x,4) + (pow(154.0,5)*par[16])/pow(x,5) + 
//         (pow(154.0,6)*par[17])/pow(x,6) + (pow(154.0,7)*par[18])/pow(x,7) + (pow(154.0,8)*par[19])/pow(x,8) + (pow(154.0,9)*par[20])/pow(x,9),2));          
//}


double coeffsecond(double *par, double x)
{
  const double r = 154.0;
  const double r2 = r*r;
  const double r3 = r*r2;
  const double r4 = r*r3;
  const double r5 = r*r4;
  const double r6 = r*r5;
  const double r7 = r*r6;
  const double r8 = r*r7;
  const double r9 = r*r8;
  const double ix = 1.0/x;
  const double ix2 = ix*ix;
  const double ix3 = ix*ix2;
  const double ix4 = ix*ix3;
  const double ix5 = ix*ix4;
  const double ix6 = ix*ix5;
  const double ix7 = ix*ix6;
  const double ix8 = ix*ix7;
  const double ix9 = ix*ix8;
  const double ix10 = ix*ix9;
  const double ix11 = ix*ix10;
   return (-2*(-((r*par[2])*ix2) - (2*r2*par[3])*ix3 - (3*r3*par[4])*ix4 - (4*r4*par[5])*ix5 - (5*r5*par[6])*ix6 - 
        (6*r6*par[7])*ix7 - (7*r7*par[8])*ix8 - (8*r8*par[9])*ix9 - (9*r9*par[10])*ix10)*
      (-((r*par[12])*ix2) - (2*r2*par[13])*ix3 - (3*r3*par[14])*ix4 - (4*r4*par[15])*ix5 - (5*r5*par[16])*ix6 - 
        (6*r6*par[17])*ix7 - (7*r7*par[18])*ix8 - (8*r8*par[19])*ix9 - (9*r9*par[20])*ix10))/
    pow(par[11] + (r*par[12])*ix + (r2*par[13])*ix2 + (r3*par[14])*ix3 + (r4*par[15])*ix4 + (r5*par[16])*ix5 + 
      (r6*par[17])*ix6 + (r7*par[18])*ix7 + (r8*par[19])*ix8 + (r9*par[20])*ix9,2) + 
   ((2*r*par[2])*ix3 + (6*r2*par[3])*ix4 + (12*r3*par[4])*ix5 + 
      (20*r4*par[5])*ix6 + (30*r5*par[6])*ix7 + (42*r6*par[7])*ix8 + (56*r7*par[8])*ix9 + (72*r8*par[9])*ix10 + 
      (90*r9*par[10])*ix11)/
    (par[11] + (r*par[12])*ix + (r2*par[13])*ix2 + (r3*par[14])*ix3 + (r4*par[15])*ix4 + (r5*par[16])*ix5 + 
      (r6*par[17])*ix6 + (r7*par[18])*ix7 + (r8*par[19])*ix8 + (r9*par[20])*ix9) + 
   (par[1] + (r*par[2])*ix + (r2*par[3])*ix2 + (r3*par[4])*ix3 + (r4*par[5])*ix4 + (r5*par[6])*ix5 + 
      (r6*par[7])*ix6 + (r7*par[8])*ix7 + (r8*par[9])*ix8 + (r9*par[10])*ix9)*
    ((2*pow(-((r*par[12])*ix2) - (2*r2*par[13])*ix3 - (3*r3*par[14])*ix4 - (4*r4*par[15])*ix5 - (5*r5*par[16])*ix6 - 
           (6*r6*par[17])*ix7 - (7*r7*par[18])*ix8 - (8*r8*par[19])*ix9 - (9*r9*par[20])*ix10,2))/
       pow(par[11] + (r*par[12])*ix + (r2*par[13])*ix2 + (r3*par[14])*ix3 + (r4*par[15])*ix4 + (r5*par[16])*ix5 + 
         (r6*par[17])*ix6 + (r7*par[18])*ix7 + (r8*par[19])*ix8 + (r9*par[20])*ix9,3) - 
      ((2*r*par[12])*ix3 + (6*r2*par[13])*ix4 + (12*r3*par[14])*ix5 + (20*r4*par[15])*ix6 + (30*r5*par[16])*ix7 + 
         (42*r6*par[17])*ix8 + (56*r7*par[18])*ix9 + (72*r8*par[19])*ix10 + (90*r9*par[20])*ix11)/
       pow(par[11] + (r*par[12])*ix + (r2*par[13])*ix2 + (r3*par[14])*ix3 + (r4*par[15])*ix4 + (r5*par[16])*ix5 + 
         (r6*par[17])*ix6 + (r7*par[18])*ix7 + (r8*par[19])*ix8 + (r9*par[20])*ix9,2));          
}


// For CHI200
double coeffMod(double *par, double x){
    return par[21] + exp(-par[1]*(200.0/x) - par[2]*(200.0/x)*(200.0/x))*par[3]*(1.0+tanh(par[4]*(x/200.0)+par[5]));
}
double coeffprimeMod(double *par, double x){
    return (1.0/x*par[1]*(200.0/x) + 2.0/x*par[2]*(200.0/x)*(200.0/x))*exp(-par[1]*(200.0/x) - par[2]*(200.0/x)*(200.0/x))*par[3]*(1.0+tanh(par[4]*(x/200.0)+par[5]))
            + exp(-par[1]*(200.0/x) - par[2]*(200.0/x)*(200.0/x))*par[3]*par[4]/200.0*pow(cosh(par[4]*(x/200.0)+par[5]),-2.0);    
}
double coeffsecondMod(double *par, double x){
    return 2.0*exp(- par[2]*(200.0/x)*(200.0/x) - par[1]*(200.0/x))*par[3]*par[4]*(2.0/x*par[2]*(200.0/x)*(200.0/x) + 1.0/x*par[1]*(200.0/x))
                  *pow(cosh(par[5] + par[4]*x/200.0),-2)/200.0 
            - 2.0*exp(- par[2]*(200.0/x)*(200.0/x) - par[1]*(200.0/x))*par[3]*par[4]*par[4]*sinh(par[5] + par[4]*x/200.0)*pow(cosh(par[5] + par[4]*x/200.0),-3)/(200.0*200.0) 
            + par[3]*exp(- par[2]*(200.0/x)*(200.0/x) - par[1]*(200.0/x))*(- 6.0/(x*x)*par[2]*(200.0/x)*(200.0/x) - 2.0/(x*x)*par[1]*(200.0/x) +pow(2.0/x*par[2]*(200.0/x)*(200.0/x) 
               + 1.0/x*par[1]*(200.0/x),2))*(1.0 + tanh(par[5] + par[4]*x/200.0));    
}


// -------------- LOW T EXTENSIONS ----------------- //
/*
double lowT_coeff(double *par, double x)
{
	return par[1]*exp(-par[3]/(x*x)-par[2]/x);
}
double lowT_coeffprime(double *par,double x)
{
	return par[1]*exp(-par[3]/(x*x)-par[2]/x)*(2.0*par[3] + par[2]*x)/(x*x*x);
}
double lowT_coeffsecond(double *par, double x)
{
	return (par[1]/(x*x*x*x*x*x))*exp(-par[3]/(x*x)-par[2]/x)
			*(4.0*par[3]*par[3] + 2.0*par[3]*(2.0*par[2] - 3.0*x)*x + par[2]*(par[2] - 2.0*x)*x*x);
}



void set_lowT_parameters(double *par1, double *par2)
{
	const double T0 = T_min_matching;
	double chi      = coeff(par1,T0);
	double chip     = coeffprime(par1,T0);
	double chipp    = coeffsecond(par1,T0);

	par2[2] = T0*T0*(3.0*chi*chip + T0*(chi*chipp - chip*chip))/(chi*chi);
	par2[3] = -(T0*T0*T0*(2.0*chi*chip + T0*(chi*chipp - chip*chip)))/(2.0*chi*chi);
	par2[1] = chi*exp( par2[3]/(T0*T0) + par2[2]/T0 );
	printf("%g %g %g %g %g\n", par2[2], par2[3], chi, chip, chipp);
}

void set_lowT_Mod_parameters(double *par1, double *par2)
{
	const double T0 = T_min_matching;
	double chi      = coeffMod(par1,T0);
	double chip     = coeffprimeMod(par1,T0);
	double chipp    = coeffsecondMod(par1,T0);

	par2[2] = T0*T0*(3.0*chi*chip + T0*(chi*chipp - chip*chip))/(chi*chi);
	par2[3] = -(T0*T0*T0*(2.0*chi*chip + T0*(chi*chipp - chip*chip)))/(2.0*chi*chi);
	par2[1] = chi*exp( par2[3]/(T0*T0) + par2[2]/T0 );
	printf("%g %g %g %g %g\n", par2[2], par2[3], chi, chip, chipp);
}
*/



double lowT_coeff(double *par, double x)
{
	return par[1]*exp(-par[2]/x);
}
double lowT_coeffprime(double *par,double x)
{
	return par[1]*par[2]*exp(-par[2]/x)/(x*x);
}
double lowT_coeffsecond(double *par, double x)
{
	return par[1]*par[2]*exp(-par[2]/x)*(par[2]-2.0*x)/(x*x*x*x);
}



void set_lowT_parameters(double *par1, double *par2)
{
	const double T0 = T_min_matching;
	double chi      = coeff(par1,T0);
	double chip     = coeffprime(par1,T0);

	par2[2] = T0*T0*chip/chi;
	par2[1] = chi*exp( par2[2]/T0 );
	//printf("%g %g %g %g\n", par2[1], par2[2], chi, chip);
}

void set_lowT_Mod_parameters(double *par1, double *par2)
{
	const double T0 = T_min_matching;
	double chi      = coeffMod(par1,T0);
	double chip     = coeffprimeMod(par1,T0);

	par2[2] = T0*T0*chip/chi;
	par2[1] = chi*exp( par2[2]/T0 );
	//printf("%g %g %g %g\n", par2[1], par2[2], chi, chip);
}



// ------------- THE COEFFICIENTS ------------------ //
double CHI000(double T){
    return coeff(CHI000PAR,T);
}
double CHI200(double T){
    return coeffMod(CHI200PAR,T);
}
double CHI020(double T){
    return coeff(CHI020PAR,T);
}
double CHI002(double T){
    return coeff(CHI002PAR,T);
}
double CHI110(double T){
    return coeff(CHI110PAR,T);
}
double CHI101(double T){
    return coeff(CHI101PAR,T);
}
double CHI011(double T){
    return coeff(CHI011PAR,T);
}
double CHI400(double T){
    return coeff(CHI400PAR,T);
}
double CHI040(double T){
    return coeff(CHI040PAR,T);
}
double CHI004(double T){
    return coeff(CHI004PAR,T);
}
double CHI310(double T){
    return coeff(CHI310PAR,T);
}
double CHI301(double T){
    return coeff(CHI301PAR,T);
}
double CHI031(double T){
    return coeff(CHI031PAR,T);
}
double CHI130(double T){
    return coeff(CHI130PAR,T);
}
double CHI103(double T){
    return coeff(CHI103PAR,T);
}
double CHI013(double T){
    return coeff(CHI013PAR,T);
}
double CHI220(double T){
    return coeff(CHI220PAR,T);
}
double CHI202(double T){
    return coeff(CHI202PAR,T);
}
double CHI022(double T){
    return coeff(CHI022PAR,T);
}
double CHI211(double T){
    return coeff(CHI211PAR,T);
}
double CHI121(double T){
    return coeff(CHI121PAR,T);
}
double CHI112(double T){
    return coeff(CHI112PAR,T);
}


// ---- Derivatives wrt T --------- //
double DCHI000DT(double T){
    return coeffprime(CHI000PAR,T);
}
double DCHI200DT(double T){
    return coeffprimeMod(CHI200PAR,T);
}
double DCHI020DT(double T){
    return coeffprime(CHI020PAR,T);
}
double DCHI002DT(double T){
    return coeffprime(CHI002PAR,T);
}
double DCHI110DT(double T){
    return coeffprime(CHI110PAR,T);
}
double DCHI101DT(double T){
    return coeffprime(CHI101PAR,T);
}
double DCHI011DT(double T){
    return coeffprime(CHI011PAR,T);
}
double DCHI400DT(double T){
    return coeffprime(CHI400PAR,T);
}
double DCHI040DT(double T){
    return coeffprime(CHI040PAR,T);
}
double DCHI004DT(double T){
    return coeffprime(CHI004PAR,T);
}
double DCHI310DT(double T){
    return coeffprime(CHI310PAR,T);
}
double DCHI301DT(double T){
    return coeffprime(CHI301PAR,T);
}
double DCHI031DT(double T){
    return coeffprime(CHI031PAR,T);
}
double DCHI130DT(double T){
    return coeffprime(CHI130PAR,T);
}
double DCHI103DT(double T){
    return coeffprime(CHI103PAR,T);
}
double DCHI013DT(double T){
    return coeffprime(CHI013PAR,T);
}
double DCHI220DT(double T){
    return coeffprime(CHI220PAR,T);
}
double DCHI202DT(double T){
    return coeffprime(CHI202PAR,T);
}
double DCHI022DT(double T){
    return coeffprime(CHI022PAR,T);
}
double DCHI211DT(double T){
    return coeffprime(CHI211PAR,T);
}
double DCHI121DT(double T){
    return coeffprime(CHI121PAR,T);
}
double DCHI112DT(double T){
    return coeffprime(CHI112PAR,T);
}


// ---- Double derivatives wrt T --------- //
double D2CHI000DT2(double T){
    return coeffsecond(CHI000PAR,T);
}
double D2CHI200DT2(double T){
    return coeffsecondMod(CHI200PAR,T);
}
double D2CHI020DT2(double T){
    return coeffsecond(CHI020PAR,T);
}
double D2CHI002DT2(double T){
    return coeffsecond(CHI002PAR,T);
}
double D2CHI110DT2(double T){
    return coeffsecond(CHI110PAR,T);
}
double D2CHI101DT2(double T){
    return coeffsecond(CHI101PAR,T);
}
double D2CHI011DT2(double T){
    return coeffsecond(CHI011PAR,T);
}
double D2CHI400DT2(double T){
    return coeffsecond(CHI400PAR,T);
}
double D2CHI040DT2(double T){
    return coeffsecond(CHI040PAR,T);
}
double D2CHI004DT2(double T){
    return coeffsecond(CHI004PAR,T);
}
double D2CHI310DT2(double T){
    return coeffsecond(CHI310PAR,T);
}
double D2CHI301DT2(double T){
    return coeffsecond(CHI301PAR,T);
}
double D2CHI031DT2(double T){
    return coeffsecond(CHI031PAR,T);
}
double D2CHI130DT2(double T){
    return coeffsecond(CHI130PAR,T);
}
double D2CHI103DT2(double T){
    return coeffsecond(CHI103PAR,T);
}
double D2CHI013DT2(double T){
    return coeffsecond(CHI013PAR,T);
}
double D2CHI220DT2(double T){
    return coeffsecond(CHI220PAR,T);
}
double D2CHI202DT2(double T){
    return coeffsecond(CHI202PAR,T);
}
double D2CHI022DT2(double T){
    return coeffsecond(CHI022PAR,T);
}
double D2CHI211DT2(double T){
    return coeffsecond(CHI211PAR,T);
}
double D2CHI121DT2(double T){
    return coeffsecond(CHI121PAR,T);
}
double D2CHI112DT2(double T){
    return coeffsecond(CHI112PAR,T);
}





///////////////////////////////////////////////////////
//           BEGIN LOW-T EXTENSION SECTION           //
///////////////////////////////////////////////////////
// ------------- THE COEFFICIENTS ------------------ //
double lowT_CHI000(double T){ return lowT_coeff(CHI000PAR_ABC,T); }
double lowT_CHI200(double T){ return lowT_coeff(CHI200PAR_ABC,T);}
double lowT_CHI020(double T){ return lowT_coeff(CHI020PAR_ABC,T); }
double lowT_CHI002(double T){ return lowT_coeff(CHI002PAR_ABC,T); }
double lowT_CHI110(double T){ return lowT_coeff(CHI110PAR_ABC,T); }
double lowT_CHI101(double T){ return lowT_coeff(CHI101PAR_ABC,T); }
double lowT_CHI011(double T){ return lowT_coeff(CHI011PAR_ABC,T); }
double lowT_CHI400(double T){ return lowT_coeff(CHI400PAR_ABC,T); }
double lowT_CHI040(double T){ return lowT_coeff(CHI040PAR_ABC,T); }
double lowT_CHI004(double T){ return lowT_coeff(CHI004PAR_ABC,T); }
double lowT_CHI310(double T){ return lowT_coeff(CHI310PAR_ABC,T); }
double lowT_CHI301(double T){ return lowT_coeff(CHI301PAR_ABC,T); }
double lowT_CHI031(double T){ return lowT_coeff(CHI031PAR_ABC,T); }
double lowT_CHI130(double T){ return lowT_coeff(CHI130PAR_ABC,T); }
double lowT_CHI103(double T){ return lowT_coeff(CHI103PAR_ABC,T); }
double lowT_CHI013(double T){ return lowT_coeff(CHI013PAR_ABC,T); }
double lowT_CHI220(double T){ return lowT_coeff(CHI220PAR_ABC,T); }
double lowT_CHI202(double T){ return lowT_coeff(CHI202PAR_ABC,T); }
double lowT_CHI022(double T){ return lowT_coeff(CHI022PAR_ABC,T); }
double lowT_CHI211(double T){ return lowT_coeff(CHI211PAR_ABC,T); }
double lowT_CHI121(double T){ return lowT_coeff(CHI121PAR_ABC,T); }
double lowT_CHI112(double T){ return lowT_coeff(CHI112PAR_ABC,T); }


// ---- Derivatives wrt T --------- //
double lowT_DCHI000DT(double T){ return lowT_coeffprime(CHI000PAR_ABC,T); }
double lowT_DCHI200DT(double T){ return lowT_coeffprime(CHI200PAR_ABC,T); }
double lowT_DCHI020DT(double T){ return lowT_coeffprime(CHI020PAR_ABC,T); }
double lowT_DCHI002DT(double T){ return lowT_coeffprime(CHI002PAR_ABC,T); }
double lowT_DCHI110DT(double T){ return lowT_coeffprime(CHI110PAR_ABC,T); }
double lowT_DCHI101DT(double T){ return lowT_coeffprime(CHI101PAR_ABC,T); }
double lowT_DCHI011DT(double T){ return lowT_coeffprime(CHI011PAR_ABC,T); }
double lowT_DCHI400DT(double T){ return lowT_coeffprime(CHI400PAR_ABC,T); }
double lowT_DCHI040DT(double T){ return lowT_coeffprime(CHI040PAR_ABC,T); }
double lowT_DCHI004DT(double T){ return lowT_coeffprime(CHI004PAR_ABC,T); }
double lowT_DCHI310DT(double T){ return lowT_coeffprime(CHI310PAR_ABC,T); }
double lowT_DCHI301DT(double T){ return lowT_coeffprime(CHI301PAR_ABC,T); }
double lowT_DCHI031DT(double T){ return lowT_coeffprime(CHI031PAR_ABC,T); }
double lowT_DCHI130DT(double T){ return lowT_coeffprime(CHI130PAR_ABC,T); }
double lowT_DCHI103DT(double T){ return lowT_coeffprime(CHI103PAR_ABC,T); }
double lowT_DCHI013DT(double T){ return lowT_coeffprime(CHI013PAR_ABC,T); }
double lowT_DCHI220DT(double T){ return lowT_coeffprime(CHI220PAR_ABC,T); }
double lowT_DCHI202DT(double T){ return lowT_coeffprime(CHI202PAR_ABC,T); }
double lowT_DCHI022DT(double T){ return lowT_coeffprime(CHI022PAR_ABC,T); }
double lowT_DCHI211DT(double T){ return lowT_coeffprime(CHI211PAR_ABC,T); }
double lowT_DCHI121DT(double T){ return lowT_coeffprime(CHI121PAR_ABC,T); }
double lowT_DCHI112DT(double T){ return lowT_coeffprime(CHI112PAR_ABC,T); }


// ---- Double derivatives wrt T --------- //
double lowT_D2CHI000DT2(double T){ return lowT_coeffsecond(CHI000PAR_ABC,T); }
double lowT_D2CHI200DT2(double T){ return lowT_coeffsecond(CHI200PAR_ABC,T); }
double lowT_D2CHI020DT2(double T){ return lowT_coeffsecond(CHI020PAR_ABC,T); }
double lowT_D2CHI002DT2(double T){ return lowT_coeffsecond(CHI002PAR_ABC,T); }
double lowT_D2CHI110DT2(double T){ return lowT_coeffsecond(CHI110PAR_ABC,T); }
double lowT_D2CHI101DT2(double T){ return lowT_coeffsecond(CHI101PAR_ABC,T); }
double lowT_D2CHI011DT2(double T){ return lowT_coeffsecond(CHI011PAR_ABC,T); }
double lowT_D2CHI400DT2(double T){ return lowT_coeffsecond(CHI400PAR_ABC,T); }
double lowT_D2CHI040DT2(double T){ return lowT_coeffsecond(CHI040PAR_ABC,T); }
double lowT_D2CHI004DT2(double T){ return lowT_coeffsecond(CHI004PAR_ABC,T); }
double lowT_D2CHI310DT2(double T){ return lowT_coeffsecond(CHI310PAR_ABC,T); }
double lowT_D2CHI301DT2(double T){ return lowT_coeffsecond(CHI301PAR_ABC,T); }
double lowT_D2CHI031DT2(double T){ return lowT_coeffsecond(CHI031PAR_ABC,T); }
double lowT_D2CHI130DT2(double T){ return lowT_coeffsecond(CHI130PAR_ABC,T); }
double lowT_D2CHI103DT2(double T){ return lowT_coeffsecond(CHI103PAR_ABC,T); }
double lowT_D2CHI013DT2(double T){ return lowT_coeffsecond(CHI013PAR_ABC,T); }
double lowT_D2CHI220DT2(double T){ return lowT_coeffsecond(CHI220PAR_ABC,T); }
double lowT_D2CHI202DT2(double T){ return lowT_coeffsecond(CHI202PAR_ABC,T); }
double lowT_D2CHI022DT2(double T){ return lowT_coeffsecond(CHI022PAR_ABC,T); }
double lowT_D2CHI211DT2(double T){ return lowT_coeffsecond(CHI211PAR_ABC,T); }
double lowT_D2CHI121DT2(double T){ return lowT_coeffsecond(CHI121PAR_ABC,T); }
double lowT_D2CHI112DT2(double T){ return lowT_coeffsecond(CHI112PAR_ABC,T); }
///////////////////////////////////////////////////////
//            END LOW-T EXTENSION SECTION            //
///////////////////////////////////////////////////////





// ------------ THERMODYNAMICS --------------------- //
double PressTaylor(double T, double muB, double muQ, double muS){
	if (T<T_min_matching)
    return lowT_CHI000(T) + lowT_CHI110(T)*muB/T*muQ/T + lowT_CHI101(T)*muB/T*muS/T
			+ lowT_CHI011(T)*muQ/T*muS/T + 1.0/2.0*lowT_CHI200(T)*muB/T*muB/T
			+ 1.0/2.0*lowT_CHI020(T)*muQ/T*muQ/T + 1.0/2.0*lowT_CHI002(T)*muS/T*muS/T 
            + 1.0/2.0*lowT_CHI211(T)*muB/T*muB/T*muQ/T*muS/T
			+ 1.0/2.0*lowT_CHI121(T)*muB/T*muQ/T*muQ/T*muS/T
			+ 1.0/2.0*lowT_CHI112(T)*muB/T*muQ/T*muS/T*muS/T
			+ 1.0/4.0*lowT_CHI220(T)*muB/T*muB/T*muQ/T*muQ/T 
            + 1.0/4.0*lowT_CHI202(T)*muB/T*muB/T*muS/T*muS/T
			+ 1.0/4.0*lowT_CHI022(T)*muQ/T*muQ/T*muS/T*muS/T
			+ 1.0/6.0*lowT_CHI310(T)*muB/T*muB/T*muB/T*muQ/T
			+ 1.0/6.0*lowT_CHI130(T)*muB/T*muQ/T*muQ/T*muQ/T 
            + 1.0/6.0*lowT_CHI301(T)*muB/T*muB/T*muB/T*muS/T
			+ 1.0/6.0*lowT_CHI103(T)*muB/T*muS/T*muS/T*muS/T
			+ 1.0/6.0*lowT_CHI031(T)*muQ/T*muQ/T*muQ/T*muS/T
			+ 1.0/6.0*lowT_CHI013(T)*muQ/T*muS/T*muS/T*muS/T 
            + 1.0/24.0*lowT_CHI400(T)*muB/T*muB/T*muB/T*muB/T
			+ 1.0/24.0*lowT_CHI040(T)*muQ/T*muQ/T*muQ/T*muQ/T
			+ 1.0/24.0*lowT_CHI004(T)*muS/T*muS/T*muS/T*muS/T;
	else
    return CHI000(T) + CHI110(T)*muB/T*muQ/T + CHI101(T)*muB/T*muS/T + CHI011(T)*muQ/T*muS/T + 1.0/2.0*CHI200(T)*muB/T*muB/T + 1.0/2.0*CHI020(T)*muQ/T*muQ/T + 1.0/2.0*CHI002(T)*muS/T*muS/T 
            + 1.0/2.0*CHI211(T)*muB/T*muB/T*muQ/T*muS/T + 1.0/2.0*CHI121(T)*muB/T*muQ/T*muQ/T*muS/T + 1.0/2.0*CHI112(T)*muB/T*muQ/T*muS/T*muS/T + 1.0/4.0*CHI220(T)*muB/T*muB/T*muQ/T*muQ/T 
            + 1.0/4.0*CHI202(T)*muB/T*muB/T*muS/T*muS/T + 1.0/4.0*CHI022(T)*muQ/T*muQ/T*muS/T*muS/T + 1.0/6.0*CHI310(T)*muB/T*muB/T*muB/T*muQ/T + 1.0/6.0*CHI130(T)*muB/T*muQ/T*muQ/T*muQ/T 
            + 1.0/6.0*CHI301(T)*muB/T*muB/T*muB/T*muS/T + 1.0/6.0*CHI103(T)*muB/T*muS/T*muS/T*muS/T + 1.0/6.0*CHI031(T)*muQ/T*muQ/T*muQ/T*muS/T + 1.0/6.0*CHI013(T)*muQ/T*muS/T*muS/T*muS/T 
            + 1.0/24.0*CHI400(T)*muB/T*muB/T*muB/T*muB/T + 1.0/24.0*CHI040(T)*muQ/T*muQ/T*muQ/T*muQ/T + 1.0/24.0*CHI004(T)*muS/T*muS/T*muS/T*muS/T;
}
double EntrTaylor(double T, double muB, double muQ, double muS){
	if (T<T_min_matching)
    return  1.0/(24.0*T*T*T)*(96.0*T*T*T*lowT_CHI000(T) + 24.0*muS*muS*T*lowT_CHI002(T) + 48.0*muQ*muS*T*lowT_CHI011(T)
			+ 24.0*muQ*muQ*T*lowT_CHI020(T) + 48.0*muB*muS*T*lowT_CHI101(T) + 48.0*muB*muQ*T*lowT_CHI110(T) 
            + 24.0*muB*muB*T*lowT_CHI200(T) + 24.0*T*T*T*T*lowT_DCHI000DT(T) + 12.0*muS*muS*T*T*lowT_DCHI002DT(T)
			+ muS*muS*muS*muS*lowT_DCHI004DT(T) + 24.0*muQ*muS*T*T*lowT_DCHI011DT(T) 
            + 4.0*muQ*muS*muS*muS*lowT_DCHI013DT(T) + 12.0*muQ*muQ*T*T*lowT_DCHI020DT(T)
			+ 6.0*muQ*muQ*muS*muS*lowT_DCHI022DT(T) + 4.0*muQ*muQ*muQ*muS*lowT_DCHI031DT(T)
			+ muQ*muQ*muQ*muQ*lowT_DCHI040DT(T) 
            + 24.0*muB*muS*T*T*lowT_DCHI101DT(T) + 4.0*muB*muS*muS*muS*lowT_DCHI103DT(T)
			+ 24.0*muB*muQ*T*T*lowT_DCHI110DT(T) + 12.0*muB*muQ*muS*muS*lowT_DCHI112DT(T)
			+ 12.0*muB*muQ*muQ*muS*lowT_DCHI121DT(T) 
            + 4.0*muB*muQ*muQ*muQ*lowT_DCHI130DT(T) + 12.0*muB*muB*T*T*lowT_DCHI200DT(T)
			+ 6.0*muB*muB*muS*muS*lowT_DCHI202DT(T) + 12.0*muB*muB*muQ*muS*lowT_DCHI211DT(T)
			+ 6.0*muB*muB*muQ*muQ*lowT_DCHI220DT(T) 
            + 4.0*muB*muB*muB*muS*lowT_DCHI301DT(T) + 4.0*muB*muB*muB*muQ*lowT_DCHI310DT(T)
			+ muB*muB*muB*muB*lowT_DCHI400DT(T));
	else
    return  1.0/(24.0*T*T*T)*(96.0*T*T*T*CHI000(T) + 24.0*muS*muS*T*CHI002(T) + 48.0*muQ*muS*T*CHI011(T)
			+ 24.0*muQ*muQ*T*CHI020(T) + 48.0*muB*muS*T*CHI101(T) + 48.0*muB*muQ*T*CHI110(T) 
            + 24.0*muB*muB*T*CHI200(T) + 24.0*T*T*T*T*DCHI000DT(T) + 12.0*muS*muS*T*T*DCHI002DT(T)
			+ muS*muS*muS*muS*DCHI004DT(T) + 24.0*muQ*muS*T*T*DCHI011DT(T) 
            + 4.0*muQ*muS*muS*muS*DCHI013DT(T) + 12.0*muQ*muQ*T*T*DCHI020DT(T)
			+ 6.0*muQ*muQ*muS*muS*DCHI022DT(T) + 4.0*muQ*muQ*muQ*muS*DCHI031DT(T)
			+ muQ*muQ*muQ*muQ*DCHI040DT(T) 
            + 24.0*muB*muS*T*T*DCHI101DT(T) + 4.0*muB*muS*muS*muS*DCHI103DT(T)
			+ 24.0*muB*muQ*T*T*DCHI110DT(T) + 12.0*muB*muQ*muS*muS*DCHI112DT(T)
			+ 12.0*muB*muQ*muQ*muS*DCHI121DT(T) 
            + 4.0*muB*muQ*muQ*muQ*DCHI130DT(T) + 12.0*muB*muB*T*T*DCHI200DT(T)
			+ 6.0*muB*muB*muS*muS*DCHI202DT(T) + 12.0*muB*muB*muQ*muS*DCHI211DT(T)
			+ 6.0*muB*muB*muQ*muQ*DCHI220DT(T) 
            + 4.0*muB*muB*muB*muS*DCHI301DT(T) + 4.0*muB*muB*muB*muQ*DCHI310DT(T)
			+ muB*muB*muB*muB*DCHI400DT(T));
}
double BarDensTaylor(double T, double muB, double muQ, double muS){
	if (T<T_min_matching)
    return lowT_CHI110(T)*muQ/T + lowT_CHI101(T)*muS/T + lowT_CHI200(T)*muB/T + lowT_CHI211(T)*muB/T*muQ/T*muS/T
			+ 1.0/2.0*lowT_CHI121(T)*muQ/T*muQ/T*muS/T + 1.0/2.0*lowT_CHI112(T)*muQ/T*muS/T*muS/T 
            + 1.0/2.0*lowT_CHI220(T)*muB/T*muQ/T*muQ/T + 1.0/2.0*lowT_CHI202(T)*muB/T*muS/T*muS/T
			+ 1.0/2.0*lowT_CHI310(T)*muB/T*muB/T*muQ/T + 1.0/6.0*lowT_CHI130(T)*muQ/T*muQ/T*muQ/T 
            + 1.0/2.0*lowT_CHI301(T)*muB/T*muB/T*muS/T + 1.0/6.0*lowT_CHI103(T)*muS/T*muS/T*muS/T
			+ 1.0/6.0*lowT_CHI400(T)*muB/T*muB/T*muB/T;
	else
    return CHI110(T)*muQ/T + CHI101(T)*muS/T + CHI200(T)*muB/T + CHI211(T)*muB/T*muQ/T*muS/T + 1.0/2.0*CHI121(T)*muQ/T*muQ/T*muS/T + 1.0/2.0*CHI112(T)*muQ/T*muS/T*muS/T 
            + 1.0/2.0*CHI220(T)*muB/T*muQ/T*muQ/T + 1.0/2.0*CHI202(T)*muB/T*muS/T*muS/T + 1.0/2.0*CHI310(T)*muB/T*muB/T*muQ/T + 1.0/6.0*CHI130(T)*muQ/T*muQ/T*muQ/T 
            + 1.0/2.0*CHI301(T)*muB/T*muB/T*muS/T + 1.0/6.0*CHI103(T)*muS/T*muS/T*muS/T + 1.0/6.0*CHI400(T)*muB/T*muB/T*muB/T;
}
double StrDensTaylor(double T, double muB, double muQ, double muS){
	if (T<T_min_matching)
    return lowT_CHI101(T)*muB/T + lowT_CHI011(T)*muQ/T + lowT_CHI002(T)*muS/T + 1.0/2.0*lowT_CHI211(T)*muB/T*muB/T*muQ/T
			+ 1.0/2.0*lowT_CHI121(T)*muB/T*muQ/T*muQ/T + lowT_CHI112(T)*muB/T*muQ/T*muS/T 
            + 1.0/2.0*lowT_CHI202(T)*muB/T*muB/T*muS/T + 1.0/2.0*lowT_CHI022(T)*muQ/T*muQ/T*muS/T
			+ 1.0/6.0*lowT_CHI301(T)*muB/T*muB/T*muB/T + 1.0/2.0*lowT_CHI103(T)*muB/T*muS/T*muS/T 
            + 1.0/6.0*lowT_CHI031(T)*muQ/T*muQ/T*muQ/T + 1.0/2.0*lowT_CHI013(T)*muQ/T*muS/T*muS/T
			+ 1.0/6.0*lowT_CHI004(T)*muS/T*muS/T*muS/T;
	else
    return CHI101(T)*muB/T + CHI011(T)*muQ/T + CHI002(T)*muS/T + 1.0/2.0*CHI211(T)*muB/T*muB/T*muQ/T
			+ 1.0/2.0*CHI121(T)*muB/T*muQ/T*muQ/T + CHI112(T)*muB/T*muQ/T*muS/T 
            + 1.0/2.0*CHI202(T)*muB/T*muB/T*muS/T + 1.0/2.0*CHI022(T)*muQ/T*muQ/T*muS/T
			+ 1.0/6.0*CHI301(T)*muB/T*muB/T*muB/T + 1.0/2.0*CHI103(T)*muB/T*muS/T*muS/T 
            + 1.0/6.0*CHI031(T)*muQ/T*muQ/T*muQ/T + 1.0/2.0*CHI013(T)*muQ/T*muS/T*muS/T
			+ 1.0/6.0*CHI004(T)*muS/T*muS/T*muS/T;
}
double ChDensTaylor(double T, double muB, double muQ, double muS){
	if (T<T_min_matching)
    return lowT_CHI110(T)*muB/T + lowT_CHI011(T)*muS/T + lowT_CHI020(T)*muQ/T + 1.0/2.0*lowT_CHI211(T)*muB/T*muB/T*muS/T
			+ lowT_CHI121(T)*muB/T*muQ/T*muS/T + 1.0/2.0*lowT_CHI112(T)*muB/T*muS/T*muS/T 
            + 1.0/2.0*lowT_CHI220(T)*muB/T*muB/T*muQ/T + 1.0/2.0*lowT_CHI022(T)*muQ/T*muS/T*muS/T
			+ 1.0/6.0*lowT_CHI310(T)*muB/T*muB/T*muB/T + 1.0/2.0*lowT_CHI130(T)*muB/T*muQ/T*muQ/T 
            + 1.0/2.0*lowT_CHI031(T)*muQ/T*muQ/T*muS/T + 1.0/6.0*lowT_CHI013(T)*muS/T*muS/T*muS/T
			+ 1.0/6.0*lowT_CHI040(T)*muQ/T*muQ/T*muQ/T;
	else
    return CHI110(T)*muB/T + CHI011(T)*muS/T + CHI020(T)*muQ/T + 1.0/2.0*CHI211(T)*muB/T*muB/T*muS/T
			+ CHI121(T)*muB/T*muQ/T*muS/T + 1.0/2.0*CHI112(T)*muB/T*muS/T*muS/T 
            + 1.0/2.0*CHI220(T)*muB/T*muB/T*muQ/T + 1.0/2.0*CHI022(T)*muQ/T*muS/T*muS/T
			+ 1.0/6.0*CHI310(T)*muB/T*muB/T*muB/T + 1.0/2.0*CHI130(T)*muB/T*muQ/T*muQ/T 
            + 1.0/2.0*CHI031(T)*muQ/T*muQ/T*muS/T + 1.0/6.0*CHI013(T)*muS/T*muS/T*muS/T
			+ 1.0/6.0*CHI040(T)*muQ/T*muQ/T*muQ/T;
}


// Higher order derivatives, Taylor expanded
//Analytical expressions
// Second order
double Chi2BTaylor(double T, double muB, double muQ, double muS){
	if (T<T_min_matching)
    return lowT_CHI200(T) + lowT_CHI211(T)*muQ/T*muS/T + 1.0/2.0*lowT_CHI220(T)*muQ/T*muQ/T
			+ 1.0/2.0*lowT_CHI202(T)*muS/T*muS/T + lowT_CHI310(T)*muB/T*muQ/T + lowT_CHI301(T)*muB/T*muS/T
			+ 1.0/2.0*lowT_CHI400(T)*muB/T*muB/T;
	else
    return CHI200(T) + CHI211(T)*muQ/T*muS/T + 1.0/2.0*CHI220(T)*muQ/T*muQ/T
			+ 1.0/2.0*CHI202(T)*muS/T*muS/T + CHI310(T)*muB/T*muQ/T + CHI301(T)*muB/T*muS/T
			+ 1.0/2.0*CHI400(T)*muB/T*muB/T;

}
double Chi2QTaylor(double T, double muB, double muQ, double muS){
	if (T<T_min_matching)
    return lowT_CHI020(T) + lowT_CHI121(T)*muB/T*muS/T + 1.0/2.0*lowT_CHI220(T)*muB/T*muB/T
			+ 1.0/2.0*lowT_CHI022(T)*muS/T*muS/T + lowT_CHI130(T)*muB/T*muQ/T + lowT_CHI031(T)*muQ/T*muS/T
			+ 1.0/2.0*lowT_CHI040(T)*muQ/T*muQ/T;
	else
    return CHI020(T) + CHI121(T)*muB/T*muS/T + 1.0/2.0*CHI220(T)*muB/T*muB/T
			+ 1.0/2.0*CHI022(T)*muS/T*muS/T + CHI130(T)*muB/T*muQ/T + CHI031(T)*muQ/T*muS/T
			+ 1.0/2.0*CHI040(T)*muQ/T*muQ/T;
}
double Chi2STaylor(double T, double muB, double muQ, double muS){
	if (T<T_min_matching)
    return lowT_CHI002(T) + lowT_CHI112(T)*muB/T*muQ/T + 1.0/2.0*lowT_CHI202(T)*muB/T*muB/T
			+ 1.0/2.0*lowT_CHI022(T)*muQ/T*muQ/T + lowT_CHI103(T)*muB/T*muS/T + lowT_CHI013(T)*muQ/T*muS/T
			+ 1.0/2.0*lowT_CHI004(T)*muS/T*muS/T;
	else
    return CHI002(T) + CHI112(T)*muB/T*muQ/T + 1.0/2.0*CHI202(T)*muB/T*muB/T
			+ 1.0/2.0*CHI022(T)*muQ/T*muQ/T + CHI103(T)*muB/T*muS/T + CHI013(T)*muQ/T*muS/T
			+ 1.0/2.0*CHI004(T)*muS/T*muS/T;
}
double Chi11BQTaylor(double T, double muB, double muQ, double muS){
	if (T<T_min_matching)
    return lowT_CHI110(T) + lowT_CHI211(T)*muB/T*muS/T + lowT_CHI121(T)*muQ/T*muS/T + 1.0/2.0*lowT_CHI112(T)*muS/T*muS/T
			+ lowT_CHI220(T)*muB/T*muQ/T + 1.0/2.0*lowT_CHI310(T)*muB/T*muB/T + 1.0/2.0*lowT_CHI130(T)*muQ/T*muQ/T;
	else
    return CHI110(T) + CHI211(T)*muB/T*muS/T + CHI121(T)*muQ/T*muS/T + 1.0/2.0*CHI112(T)*muS/T*muS/T
			+ CHI220(T)*muB/T*muQ/T + 1.0/2.0*CHI310(T)*muB/T*muB/T + 1.0/2.0*CHI130(T)*muQ/T*muQ/T;
}
double Chi11BSTaylor(double T, double muB, double muQ, double muS){
	if (T<T_min_matching)
    return lowT_CHI101(T) + lowT_CHI211(T)*muB/T*muQ/T + 1.0/2.0*lowT_CHI121(T)*muQ/T*muQ/T + lowT_CHI112(T)*muQ/T*muS/T
			+ lowT_CHI202(T)*muB/T*muS/T + 1.0/2.0*lowT_CHI301(T)*muB/T*muB/T + 1.0/2.0*lowT_CHI103(T)*muS/T*muS/T;
	else
    return CHI101(T) + CHI211(T)*muB/T*muQ/T + 1.0/2.0*CHI121(T)*muQ/T*muQ/T + CHI112(T)*muQ/T*muS/T
			+ CHI202(T)*muB/T*muS/T + 1.0/2.0*CHI301(T)*muB/T*muB/T + 1.0/2.0*CHI103(T)*muS/T*muS/T;
}
double Chi11QSTaylor(double T, double muB, double muQ, double muS){
	if (T<T_min_matching)
    return lowT_CHI011(T) + 1.0/2.0*lowT_CHI211(T)*muB/T*muB/T + lowT_CHI121(T)*muB/T*muQ/T + lowT_CHI112(T)*muB/T*muS/T
			+ lowT_CHI022(T)*muQ/T*muS/T + 1.0/2.0*lowT_CHI031(T)*muQ/T*muQ/T + 1.0/2.0*lowT_CHI013(T)*muS/T*muS/T;
	else
    return CHI011(T) + 1.0/2.0*CHI211(T)*muB/T*muB/T + CHI121(T)*muB/T*muQ/T + CHI112(T)*muB/T*muS/T
			+ CHI022(T)*muQ/T*muS/T + 1.0/2.0*CHI031(T)*muQ/T*muQ/T + 1.0/2.0*CHI013(T)*muS/T*muS/T;
}

// Second order (one is T)
double DBarDensDTTaylor(double T, double muB, double muQ, double muS){
	if (T<T_min_matching)
    return 1.0/(6.0*T*T)*(12.0*muS*T*lowT_CHI101(T) + 12.0*muQ*T*lowT_CHI110(T) + 12.0*muB*T*lowT_CHI200(T)
			+ 6.0*muS*T*T*lowT_DCHI101DT(T) + muS*muS*muS*lowT_DCHI103DT(T) + 6.0*muQ*T*T*lowT_DCHI110DT(T) 
            + 3.0*muQ*muS*muS*lowT_DCHI112DT(T) + 3.0*muQ*muQ*muS*lowT_DCHI121DT(T) + muQ*muQ*muQ*lowT_DCHI130DT(T)
			+ 6.0*muB*T*T*lowT_DCHI200DT(T) + 3.0*muB*muS*muS*lowT_DCHI202DT(T) 
            + 6.0*muB*muQ*muS*lowT_DCHI211DT(T) + 3.0*muB*muQ*muQ*lowT_DCHI220DT(T) + 3.0*muB*muB*muS*lowT_DCHI301DT(T)
			+ 3.0*muB*muB*muQ*lowT_DCHI310DT(T) + muB*muB*muB*lowT_DCHI400DT(T));
	else
    return 1.0/(6.0*T*T)*(12.0*muS*T*CHI101(T) + 12.0*muQ*T*CHI110(T) + 12.0*muB*T*CHI200(T)
			+ 6.0*muS*T*T*DCHI101DT(T) + muS*muS*muS*DCHI103DT(T) + 6.0*muQ*T*T*DCHI110DT(T) 
            + 3.0*muQ*muS*muS*DCHI112DT(T) + 3.0*muQ*muQ*muS*DCHI121DT(T) + muQ*muQ*muQ*DCHI130DT(T)
			+ 6.0*muB*T*T*DCHI200DT(T) + 3.0*muB*muS*muS*DCHI202DT(T) 
            + 6.0*muB*muQ*muS*DCHI211DT(T) + 3.0*muB*muQ*muQ*DCHI220DT(T) + 3.0*muB*muB*muS*DCHI301DT(T)
			+ 3.0*muB*muB*muQ*DCHI310DT(T) + muB*muB*muB*DCHI400DT(T));
}
double DStrDensDTTaylor(double T, double muB, double muQ, double muS){
	if (T<T_min_matching)
    return 1.0/(6.0*T*T)*(12.0*muS*T*lowT_CHI002(T) + 12.0*muQ*T*lowT_CHI011(T) + 12.0*muB*T*lowT_CHI101(T)
			+ 6.0*muS*T*T*lowT_DCHI002DT(T) + muS*muS*muS*lowT_DCHI004DT(T) + 6.0*muQ*T*T*lowT_DCHI011DT(T) 
            + 3.0*muQ*muS*muS*lowT_DCHI013DT(T) + 3.0*muQ*muQ*muS*lowT_DCHI022DT(T) + muQ*muQ*muQ*lowT_DCHI031DT(T)
			+ 6.0*muB*T*T*lowT_DCHI101DT(T) + 3.0*muB*muS*muS*lowT_DCHI103DT(T) 
            + 6.0*muB*muQ*muS*lowT_DCHI112DT(T) + 3.0*muB*muQ*muQ*lowT_DCHI121DT(T) + 3.0*muB*muB*muS*lowT_DCHI202DT(T)
			+ 3.0*muB*muB*muQ*lowT_DCHI211DT(T) + muB*muB*muB*lowT_DCHI301DT(T));
	else
    return 1.0/(6.0*T*T)*(12.0*muS*T*CHI002(T) + 12.0*muQ*T*CHI011(T) + 12.0*muB*T*CHI101(T)
			+ 6.0*muS*T*T*DCHI002DT(T) + muS*muS*muS*DCHI004DT(T) + 6.0*muQ*T*T*DCHI011DT(T) 
            + 3.0*muQ*muS*muS*DCHI013DT(T) + 3.0*muQ*muQ*muS*DCHI022DT(T) + muQ*muQ*muQ*DCHI031DT(T)
			+ 6.0*muB*T*T*DCHI101DT(T) + 3.0*muB*muS*muS*DCHI103DT(T) 
            + 6.0*muB*muQ*muS*DCHI112DT(T) + 3.0*muB*muQ*muQ*DCHI121DT(T) + 3.0*muB*muB*muS*DCHI202DT(T)
			+ 3.0*muB*muB*muQ*DCHI211DT(T) + muB*muB*muB*DCHI301DT(T));
}
double DChDensDTTaylor(double T, double muB, double muQ, double muS){
	if (T<T_min_matching)
    return 1.0/(6.0*T*T)*(12.0*muS*T*lowT_CHI011(T) + 12.0*muQ*T*lowT_CHI020(T) + 12.0*muB*T*lowT_CHI110(T)
			+ 6.0*muS*T*T*lowT_DCHI011DT(T) + muS*muS*muS*lowT_DCHI013DT(T) + 6.0*muQ*T*T*lowT_DCHI020DT(T) 
            + 3.0*muQ*muS*muS*lowT_DCHI022DT(T) + 3.0*muQ*muQ*muS*lowT_DCHI031DT(T) + muQ*muQ*muQ*lowT_DCHI040DT(T)
			+ 6.0*muB*T*T*lowT_DCHI110DT(T) + 3.0*muB*muS*muS*lowT_DCHI112DT(T) 
            + 6.0*muB*muQ*muS*lowT_DCHI121DT(T) + 3.0*muB*muQ*muQ*lowT_DCHI130DT(T) + 3.0*muB*muB*muS*lowT_DCHI211DT(T)
			+ 3.0*muB*muB*muQ*lowT_DCHI220DT(T) + muB*muB*muB*lowT_DCHI310DT(T));
	else
    return 1.0/(6.0*T*T)*(12.0*muS*T*CHI011(T) + 12.0*muQ*T*CHI020(T) + 12.0*muB*T*CHI110(T)
			+ 6.0*muS*T*T*DCHI011DT(T) + muS*muS*muS*DCHI013DT(T) + 6.0*muQ*T*T*DCHI020DT(T) 
            + 3.0*muQ*muS*muS*DCHI022DT(T) + 3.0*muQ*muQ*muS*DCHI031DT(T) + muQ*muQ*muQ*DCHI040DT(T)
			+ 6.0*muB*T*T*DCHI110DT(T) + 3.0*muB*muS*muS*DCHI112DT(T) 
            + 6.0*muB*muQ*muS*DCHI121DT(T) + 3.0*muB*muQ*muQ*DCHI130DT(T) + 3.0*muB*muB*muS*DCHI211DT(T)
			+ 3.0*muB*muB*muQ*DCHI220DT(T) + muB*muB*muB*DCHI310DT(T));
}

// Second order (both are T)
double DEntrDTTaylor(double T, double muB, double muQ, double muS){
	if (T<T_min_matching)
    return (288*pow(T,2)*lowT_CHI000(T) + 24*(pow(muS,2)*lowT_CHI002(T)
				+ 2*muQ*muS*lowT_CHI011(T) + pow(muQ,2)*lowT_CHI020(T)
				+ 2*muB*muS*lowT_CHI101(T) + 2*muB*muQ*lowT_CHI110(T) + pow(muB,2)*lowT_CHI200(T)) 
				+ 192*pow(T,3)*lowT_DCHI000DT(T) + 48*pow(muS,2)*T*lowT_DCHI002DT(T) + 96*muQ*muS*T*lowT_DCHI011DT(T)
				+ 48*pow(muQ,2)*T*lowT_DCHI020DT(T) + 96*muB*muS*T*lowT_DCHI101DT(T) + 96*muB*muQ*T*lowT_DCHI110DT(T) 
				+ 48*pow(muB,2)*T*lowT_DCHI200DT(T) + 24*pow(T,4)*lowT_D2CHI000DT2(T)
				+ 12*pow(muS,2)*pow(T,2)*lowT_D2CHI002DT2(T) + pow(muS,4)*lowT_D2CHI004DT2(T)
				+ 24*muQ*muS*pow(T,2)*lowT_D2CHI011DT2(T) 
				+ 4*muQ*pow(muS,3)*lowT_D2CHI013DT2(T) + 12*pow(muQ,2)*pow(T,2)*lowT_D2CHI020DT2(T)
				+ 6*pow(muQ,2)*pow(muS,2)*lowT_D2CHI022DT2(T) + 4*pow(muQ,3)*muS*lowT_D2CHI031DT2(T)
				+ pow(muQ,4)*lowT_D2CHI040DT2(T) 
				+ 24*muB*muS*pow(T,2)*lowT_D2CHI101DT2(T) + 4*muB*pow(muS,3)*lowT_D2CHI103DT2(T)
				+ 24*muB*muQ*pow(T,2)*lowT_D2CHI110DT2(T) + 12*muB*muQ*pow(muS,2)*lowT_D2CHI112DT2(T) 
				+ 12*muB*pow(muQ,2)*muS*lowT_D2CHI121DT2(T) + 4*muB*pow(muQ,3)*lowT_D2CHI130DT2(T)
				+ 12*pow(muB,2)*pow(T,2)*lowT_D2CHI200DT2(T) + 6*pow(muB,2)*pow(muS,2)*lowT_D2CHI202DT2(T) 
				+ 12*pow(muB,2)*muQ*muS*lowT_D2CHI211DT2(T) + 6*pow(muB,2)*pow(muQ,2)*lowT_D2CHI220DT2(T)
				+ 4*pow(muB,3)*muS*lowT_D2CHI301DT2(T) + 4*pow(muB,3)*muQ*lowT_D2CHI310DT2(T) 
				+ pow(muB,4)*lowT_D2CHI400DT2(T))/(24.*pow(T,2));
	else
    return (288*pow(T,2)*CHI000(T) + 24*(pow(muS,2)*CHI002(T) + 2*muQ*muS*CHI011(T) + pow(muQ,2)*CHI020(T) + 2*muB*muS*CHI101(T) + 2*muB*muQ*CHI110(T) + pow(muB,2)*CHI200(T)) 
               + 192*pow(T,3)*DCHI000DT(T) + 48*pow(muS,2)*T*DCHI002DT(T) + 96*muQ*muS*T*DCHI011DT(T) + 48*pow(muQ,2)*T*DCHI020DT(T) + 96*muB*muS*T*DCHI101DT(T) + 96*muB*muQ*T*DCHI110DT(T) 
               + 48*pow(muB,2)*T*DCHI200DT(T) + 24*pow(T,4)*D2CHI000DT2(T) + 12*pow(muS,2)*pow(T,2)*D2CHI002DT2(T) + pow(muS,4)*D2CHI004DT2(T) + 24*muQ*muS*pow(T,2)*D2CHI011DT2(T) 
               + 4*muQ*pow(muS,3)*D2CHI013DT2(T) + 12*pow(muQ,2)*pow(T,2)*D2CHI020DT2(T) + 6*pow(muQ,2)*pow(muS,2)*D2CHI022DT2(T) + 4*pow(muQ,3)*muS*D2CHI031DT2(T) + pow(muQ,4)*D2CHI040DT2(T) 
               + 24*muB*muS*pow(T,2)*D2CHI101DT2(T) + 4*muB*pow(muS,3)*D2CHI103DT2(T) + 24*muB*muQ*pow(T,2)*D2CHI110DT2(T) + 12*muB*muQ*pow(muS,2)*D2CHI112DT2(T) 
               + 12*muB*pow(muQ,2)*muS*D2CHI121DT2(T) + 4*muB*pow(muQ,3)*D2CHI130DT2(T) + 12*pow(muB,2)*pow(T,2)*D2CHI200DT2(T) + 6*pow(muB,2)*pow(muS,2)*D2CHI202DT2(T) 
               + 12*pow(muB,2)*muQ*muS*D2CHI211DT2(T) + 6*pow(muB,2)*pow(muQ,2)*D2CHI220DT2(T) + 4*pow(muB,3)*muS*D2CHI301DT2(T) + 4*pow(muB,3)*muQ*D2CHI310DT2(T) 
               + pow(muB,4)*D2CHI400DT2(T))/(24.*pow(T,2));
}







// Speed of Sound expression.
double SpSound(double T, double muB, double muQ, double muS){
   C1B = BarDensTaylor(T,muB,muS,muQ);
   C1Q = ChDensTaylor(T,muB,muS,muQ);
   C1S = StrDensTaylor(T,muB,muS,muQ);
   C1T = EntrTaylor(T,muB,muS,muQ);
   
   C2B2 = Chi2BTaylor(T,muB,muS,muQ);   
   C2Q2 = Chi2QTaylor(T,muB,muS,muQ);
   C2S2 = Chi2STaylor(T,muB,muS,muQ);
   C2BQ = Chi11BQTaylor(T,muB,muS,muQ);
   C2BS = Chi11BSTaylor(T,muB,muS,muQ);
   C2QS = Chi11QSTaylor(T,muB,muS,muQ);
   C2TB = DBarDensDTTaylor(T,muB,muS,muQ);
   C2TQ = DChDensDTTaylor(T,muB,muS,muQ);
   C2TS = DStrDensDTTaylor(T,muB,muS,muQ);
   C2T2 = DEntrDTTaylor(T,muB,muS,muQ);

//   C1B = BarDensTaylor(T,muB,muQ,muS);
//   C1Q = ChDensTaylor(T,muB,muQ,muS);
//   C1S = StrDensTaylor(T,muB,muQ,muS);
//   C1T = EntrTaylor(T,muB,muQ,muS);
   
//   C2B2 = Chi2BTaylor(T,muB,muQ,muS);   
//   C2Q2 = Chi2QTaylor(T,muB,muQ,muS);
//   C2S2 = Chi2STaylor(T,muB,muQ,muS);
//   C2BQ = Chi11BQTaylor(T,muB,muQ,muS);
//   C2BS = Chi11BSTaylor(T,muB,muQ,muS);
//   C2QS = Chi11QSTaylor(T,muB,muQ,muS);
//   C2TB = DBarDensDTTaylor(T,muB,muQ,muS);
//   C2TQ = DChDensDTTaylor(T,muB,muQ,muS);
//   C2TS = DStrDensDTTaylor(T,muB,muQ,muS);
//   C2T2 = DEntrDTTaylor(T,muB,muQ,muS);

//  printf("Check input(C): %lf  %lf  %lf  %lf  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f  %3.12f\n",
//          T, muB, muQ, muS, C1T, C1B, C1S, C1Q, C2B2, C2Q2, C2S2, C2BQ, C2BS, C2QS, C2TB, C2TQ, C2TS, C2T2);
      
   return T*(-(C2BQ*C2S2*C2TQ*C1B) - C2BQ*C2S2*C2TB*C1Q - pow(C2BS,2)*C2TQ*C1Q + C2B2*C2S2*C2TQ*C1Q + C2BQ*C2BS*C2TS*C1Q + C2BQ*C2BS*C2TQ*C1S - pow(C2BQ,2)*C2TS*C1S + pow(C2BQ,2)*C2S2*C1T 
            + pow(C2QS,2)*(-(C2TB*C1B) + C2B2*C1T) + C2Q2*(C2S2*C2TB*C1B - C2BS*C2TS*C1B - C2BS*C2TB*C1S + C2B2*C2TS*C1S + pow(C2BS,2)*C1T - C2B2*C2S2*C1T) + C2QS*(C2BQ*C2TS*C1B - C2B2*C2TS*C1Q 
            + C2BQ*C2TB*C1S - C2B2*C2TQ*C1S + C2BS*(C2TQ*C1B + C2TB*C1Q - 2.0*C2BQ*C1T)))
               *1.0/((pow(C2BQ,2)*C2S2*C2T2 - pow(C2QS,2)*pow(C2TB,2) + C2Q2*C2S2*pow(C2TB,2) - 2.0*C2BQ*C2S2*C2TB*C2TQ + pow(C2BS,2)*(C2Q2*C2T2 - pow(C2TQ,2)) 
            + 2.0*C2BQ*C2QS*C2TB*C2TS - pow(C2BQ,2)*pow(C2TS,2) + 2.0*C2BS*(-(C2BQ*C2QS*C2T2) + C2QS*C2TB*C2TQ - C2Q2*C2TB*C2TS + C2BQ*C2TQ*C2TS) + C2B2*(pow(C2QS,2)*C2T2 - C2Q2*C2S2*C2T2 
            + C2S2*pow(C2TQ,2) - 2.0*C2QS*C2TQ*C2TS + C2Q2*pow(C2TS,2)))*T) 
          + T*(C1B/(muB*C1B + muQ*C1Q + muS*C1S + T*C1T))*(-(C2QS*muB*(C2BQ*(C2TS*C1B + C2TB*C1S) - C2B2*(C2TS*C1Q + C2TQ*C1S) + C2BS*(C2TQ*C1B + C2TB*C1Q - 2*C2BQ*C1T))) + muB*((pow(C2BS,2) 
            - C2B2*C2S2)*C2TQ*C1Q + C2BQ*(C2S2*C2TQ*C1B + C2S2*C2TB*C1Q - C2BS*C2TS*C1Q - C2BS*C2TQ*C1S) + pow(C2BQ,2)*(C2TS*C1S - C2S2*C1T) + C2Q2*(-(C2S2*C2TB*C1B) + C2BS*C2TS*C1B 
            + C2BS*C2TB*C1S - C2B2*C2TS*C1S - pow(C2BS,2)*C1T + C2B2*C2S2*C1T)) + C2QS*(-2*C2TQ*C2TS*C1B + C2TB*C2TS*C1Q + C2TB*C2TQ*C1S - C2T2*(C2BS*C1Q + C2BQ*C1S) + C2BS*C2TQ*C1T 
            + C2BQ*C2TS*C1T)*T + (-((C2BS*C2TQ - C2BQ*C2TS)*(-(C2TS*C1Q) + C2TQ*C1S)) + C2S2*(C2BQ*C2T2*C1Q + C2TQ*(C2TQ*C1B - C2TB*C1Q - C2BQ*C1T)) + C2Q2*(C2BS*C2T2*C1S + C2TS*(C2TS*C1B 
            - C2TB*C1S - C2BS*C1T) + C2S2*(-(C2T2*C1B) + C2TB*C1T)))*T + pow(C2QS,2)*(C2TB*muB*C1B - C2B2*muB*C1T + C2T2*C1B*T - C2TB*C1T*T))/((pow(C2BQ,2)*C2S2*C2T2 - pow(C2QS,2)*pow(C2TB,2) 
            + C2Q2*C2S2*pow(C2TB,2) - 2*C2BQ*C2S2*C2TB*C2TQ + pow(C2BS,2)*(C2Q2*C2T2 - pow(C2TQ,2)) + 2*C2BQ*C2QS*C2TB*C2TS - pow(C2BQ,2)*pow(C2TS,2) + 2*C2BS*(-(C2BQ*C2QS*C2T2) + C2QS*C2TB*C2TQ 
            - C2Q2*C2TB*C2TS + C2BQ*C2TQ*C2TS) + C2B2*(pow(C2QS,2)*C2T2 - C2Q2*C2S2*C2T2 + C2S2*pow(C2TQ,2) - 2*C2QS*C2TQ*C2TS + C2Q2*pow(C2TS,2)))*T)
          + T*(C1Q/(muB*C1B + muQ*C1Q + muS*C1S + T*C1T))*(pow(C2QS,2)*muQ*(C2TB*C1B - C2B2*C1T) - C2QS*muQ*(C2BQ*(C2TS*C1B + C2TB*C1S) - C2B2*(C2TS*C1Q + C2TQ*C1S) + C2BS*(C2TQ*C1B + C2TB*C1Q 
            - 2*C2BQ*C1T)) + muQ*((pow(C2BS,2) - C2B2*C2S2)*C2TQ*C1Q + C2BQ*(C2S2*C2TQ*C1B + C2S2*C2TB*C1Q - C2BS*C2TS*C1Q - C2BS*C2TQ*C1S) + pow(C2BQ,2)*(C2TS*C1S - C2S2*C1T) + C2Q2*(-(C2S2*C2TB*C1B) 
            + C2BS*C2TS*C1B + C2BS*C2TB*C1S - C2B2*C2TS*C1S - pow(C2BS,2)*C1T + C2B2*C2S2*C1T)) + C2QS*(-(C2BS*C2T2*C1B) + C2TB*C2TS*C1B + C2B2*C2T2*C1S - pow(C2TB,2)*C1S + C2BS*C2TB*C1T 
            - C2B2*C2TS*C1T)*T + (C2BS*C2TQ*C2TS*C1B + pow(C2BS,2)*C2T2*C1Q - 2*C2BS*C2TB*C2TS*C1Q + C2B2*pow(C2TS,2)*C1Q + C2BS*C2TB*C2TQ*C1S - C2B2*C2TQ*C2TS*C1S - pow(C2BS,2)*C2TQ*C1T 
            + C2S2*(-(C2TB*C2TQ*C1B) - C2B2*C2T2*C1Q + pow(C2TB,2)*C1Q + C2B2*C2TQ*C1T) + C2BQ*(-(C2BS*C2T2*C1S) + C2TS*(-(C2TS*C1B) + C2TB*C1S + C2BS*C1T) + C2S2*(C2T2*C1B - C2TB*C1T)))*T)
               *1.0/((pow(C2BQ,2)*C2S2*C2T2 - pow(C2QS,2)*pow(C2TB,2) + C2Q2*C2S2*pow(C2TB,2) - 2*C2BQ*C2S2*C2TB*C2TQ + pow(C2BS,2)*(C2Q2*C2T2 - pow(C2TQ,2)) + 2*C2BQ*C2QS*C2TB*C2TS 
            - pow(C2BQ,2)*pow(C2TS,2) + 2*C2BS*(-(C2BQ*C2QS*C2T2) + C2QS*C2TB*C2TQ - C2Q2*C2TB*C2TS + C2BQ*C2TQ*C2TS) + C2B2*(pow(C2QS,2)*C2T2 - C2Q2*C2S2*C2T2 + C2S2*pow(C2TQ,2) 
            - 2*C2QS*C2TQ*C2TS + C2Q2*pow(C2TS,2)))*T)
          + T*(C1S/(muB*C1B + muQ*C1Q + muS*C1S + T*C1T))*(pow(C2QS,2)*muS*(C2TB*C1B - C2B2*C1T) - C2QS*muS*(C2BQ*(C2TS*C1B + C2TB*C1S) - C2B2*(C2TS*C1Q + C2TQ*C1S) + C2BS*(C2TQ*C1B + C2TB*C1Q 
            - 2*C2BQ*C1T)) + muS*((pow(C2BS,2) - C2B2*C2S2)*C2TQ*C1Q + C2BQ*(C2S2*C2TQ*C1B + C2S2*C2TB*C1Q - C2BS*C2TS*C1Q - C2BS*C2TQ*C1S) + pow(C2BQ,2)*(C2TS*C1S - C2S2*C1T) 
            + C2Q2*(-(C2S2*C2TB*C1B) + C2BS*C2TS*C1B + C2BS*C2TB*C1S - C2B2*C2TS*C1S - pow(C2BS,2)*C1T + C2B2*C2S2*C1T)) + C2QS*(-(C2BQ*C2T2*C1B) + C2TB*C2TQ*C1B + C2B2*C2T2*C1Q - pow(C2TB,2)*C1Q 
            + C2BQ*C2TB*C1T - C2B2*C2TQ*C1T)*T + (C2BQ*C2TQ*C2TS*C1B + C2BQ*C2TB*C2TS*C1Q - C2B2*C2TQ*C2TS*C1Q + pow(C2BQ,2)*C2T2*C1S - 2*C2BQ*C2TB*C2TQ*C1S + C2B2*pow(C2TQ,2)*C1S 
            - pow(C2BQ,2)*C2TS*C1T + C2Q2*(-(C2TB*C2TS*C1B) - C2B2*C2T2*C1S + pow(C2TB,2)*C1S + C2B2*C2TS*C1T) + C2BS*(-(C2BQ*C2T2*C1Q) + C2TQ*(-(C2TQ*C1B) + C2TB*C1Q + C2BQ*C1T) 
            + C2Q2*(C2T2*C1B - C2TB*C1T)))*T)/((pow(C2BQ,2)*C2S2*C2T2 - pow(C2QS,2)*pow(C2TB,2) + C2Q2*C2S2*pow(C2TB,2) - 2*C2BQ*C2S2*C2TB*C2TQ + pow(C2BS,2)*(C2Q2*C2T2 - pow(C2TQ,2)) 
            + 2*C2BQ*C2QS*C2TB*C2TS - pow(C2BQ,2)*pow(C2TS,2) + 2*C2BS*(-(C2BQ*C2QS*C2T2) + C2QS*C2TB*C2TQ - C2Q2*C2TB*C2TS + C2BQ*C2TQ*C2TS) + C2B2*(pow(C2QS,2)*C2T2 - C2Q2*C2S2*C2T2 
            + C2S2*pow(C2TQ,2) - 2*C2QS*C2TQ*C2TS + C2Q2*pow(C2TS,2)))*T);
}

double SpSoundReadable(double T, double muB, double muQ, double muS)
{
  C1B = BarDensTaylor(T,muB,muS,muQ);
  C1Q = ChDensTaylor(T,muB,muS,muQ);
  C1S = StrDensTaylor(T,muB,muS,muQ);
  C1T = EntrTaylor(T,muB,muS,muQ);

  double C2BB = Chi2BTaylor(T,muB,muS,muQ);   
  double C2QQ = Chi2QTaylor(T,muB,muS,muQ);
  double C2SS = Chi2STaylor(T,muB,muS,muQ);
  C2BQ = Chi11BQTaylor(T,muB,muS,muQ);
  C2BS = Chi11BSTaylor(T,muB,muS,muQ);
  C2QS = Chi11QSTaylor(T,muB,muS,muQ);
  C2TB = DBarDensDTTaylor(T,muB,muS,muQ);
  C2TQ = DChDensDTTaylor(T,muB,muS,muQ);
  C2TS = DStrDensDTTaylor(T,muB,muS,muQ);
  C2T2 = DEntrDTTaylor(T,muB,muS,muQ);

  double Delta = C2BS*C2BS*C2QQ + C2SQ*C2SQ*C2BB + C2BQ*C2BQ*C2SS
                  - 2.0*C2BQ*C2BS*C2SQ - C2BB*C2SS*C2QQ;

  double DTSSQ = C2SQ*C2TS - C2SS*C2TQ;
  double DTQQS = C2SQ*C2TQ - C2QQ*C2TS;
  double DSQSQ = C2QQ*C2SS - C2SQ*C2SQ;
  double DTQQB = C2BQ*C2TQ - C2QQ*C2TB;
  double DTBBQ = C2BQ*C2TB - C2BB*C2TQ;
  double DBQBQ = C2BB*C2QQ - C2BQ*C2BQ;
  double DTSSB = C2BS*C2TS - C2SS*C2TB;
  double DTBBS = C2BS*C2TB - C2BB*C2TS;
  double DBSBS = C2BB*C2SS - C2BS*C2BS;
  double DTBQS = C2BS*C2TQ - C2BQ*C2TS;
  double DTSQB = C2BS*C2TQ - C2SQ*C2TB;
  double DBQQS = C2BQ*C2SQ - C2BS*C2QQ;
  double DTQQS = C2SQ*C2TQ - C2QQ*C2TS;
  double DTQQB = C2BQ*C2TQ - C2QQ*C2TB;
  double DTBSQ = C2BQ*C2TS - C2BS*C2TQ;
  double DTQSB = C2BQ*C2TS - C2SQ*C2TB;
  double DBSSQ = C2BS*C2SQ - C2BQ*C2SS;
  double DTSSQ = C2SQ*C2TS - C2SS*C2TQ;
  double DTSSB = C2BS*C2TS - C2SS*C2TB;
  double DTQBS = C2SQ*C2TB - C2BQ*C2TS;
  double DTSBQ = C2SQ*C2TB - C2BS*C2TQ;
  double DSBBQ = C2BQ*C2BS - C2BB*C2SQ;
  double DTBBQ = C2BQ*C2TB - C2BB*C2TQ;
  double DTBBS = C2BS*C2TB - C2BB*C2TS;
  double DSQSQ = C2QQ*C2SS - C2SQ*C2SQ;
  double DBSSQ = C2BS*C2SQ - C2BQ*C2SS;
  double DBQQS = C2BQ*C2SQ - C2BS*C2QQ;
  double DBQBQ = C2BB*C2QQ - C2BQ*C2BQ;
  double DBQQS = C2BQ*C2SQ - C2BS*C2QQ;
  double DBSQB = C2BQ*C2BS - C2BB*C2SQ;
  double DBSBS = C2BB*C2SS - C2BS*C2BS;
  double DBSSQ = C2BS*C2SQ - C2BQ*C2SS;
  double DBSQB = C2BQ*C2BS - C2BB*C2SQ;

  double D1 = C2TS*DBQQS + C2TQ*DBSSQ + C2TB*DSQSQ;
  double D2 = C2TS*DBQBQ + C2TB*DBQQS + C2TQ*DBSQB;
  double D3 = C2TQ*DBSBS + C2TS*DBSQB + C2TB*DBSSQ;

  double DBB = C2TT*DSQSQ + C2TS*DTQQS + C2TQ*DTSSQ;
  double DQQ = C2TT*DBQBQ + C2TQ*DTBBQ + C2TB*DTQQB;
  double DSS = C2TT*DBSBS + C2TS*DTBBS + C2TB*DTSSB;
  double DBS = C2TT*DBQQS - 0.5*(C2TS*DTQQB + C2TB*DTQQS) + 0.5*(C2TQ*DTBQS + DTSQB);
  double DBQ = C2TT*DBSSQ - 0.5*(C2TQ*DTSSB + C2TB*DTSSQ) + 0.5*(C2TS*DTBSQ + DTQSB);
  double DSQ = C2TT*DSBBQ - 0.5*(C2TS*DTBBQ + C2TQ*DTBBS) + 0.5*(C2TB*DTQBS + DTSBQ);

  double DT = Delta;
  double DB = Delta * muB - T*D1;
  double DS = Delta * muS - T*D2;
  double DQ = Delta * muQ - T*D3;

  double MB = D1;
  double MS = D2;
  double MQ = D3;
  double MBB = T*DBB + muB*D1;
  double MBS = T*DBS + muB*D2;
  double MBQ = T*DBQ + muB*D3;
  double MSB = T*DBS + muS*D1;
  double MSS = T*DSS + muS*D2;
  double MSQ = T*DSQ + muS*D3;
  double MQB = T*DBQ + muQ*D1;
  double MQS = T*DSQ + muQ*D2;
  double MQQ = T*DQQ + muQ*D3;

  // Term 1 (dp/de at constant density).
  double dedT = T*C2TT + muB*C2TB + muS*C2TS + muQ*C2TQ;
  double dedmuB = T*C2TB + muB*C2BB + muS*C2BS + muQ*C2BQ;
  double dedmuS = T*C2TS + muB*C2BS + muS*C2SS + muQ*C2SQ;
  double dedmuQ = T*C2TQ + muB*C2BQ + muS*C2SQ + muQ*C2QQ;
  double dpde = (C1T*DT + C1B*MB + C1S*MS + C1Q*MQ)
                / (dedT*DT + dedmuB*MB + dedmuS*MS + dedmuQ*MQ);

  // Term 2 (dp/dB at constant e, S, Q).
  double dpdB = (C1T*DB + C1B*MBB + C1S*MBS + C1Q*MBQ)
                / (C2TB*DB + C2BB*MBB + C2BS*MBS + C2BQ*MBQ);

  // Term 3 (dp/dS at constant e, B, Q).
  double dpdS = (C1T*DS + C1B*MBS + C1S*MSS + C1Q*MSQ)
                / (C2TS*DS + C2BS*MBS + C2SS*MSS + C2SQ*MSQ);

  // Term 4 (dp/dQ at constant e, B, S).
  double dpdQ = (C1T*DQ + C1B*MBQ + C1S*MSQ + C1Q*MQQ)
                / (C2TQ*DQ + C2BQ*MBQ + C2SQ*MSQ + C2QQ*MQQ);

  double e_plus_p = C1T + (C1B*muB + C1S*muS + C1Q*muQ)/T;

  return dpde + (C1B*dpdB + C1S*dpdS + C1Q*dpdQ) / e_plus_p;
}


//Output derivatives

double P2B2(double T, double muB, double muQ, double muS)
{   
   return Chi2BTaylor(T,muB,muS,muQ);
}

double P2Q2(double T, double muB, double muQ, double muS){
   return Chi2QTaylor(T,muB,muS,muQ);
}

double P2S2(double T, double muB, double muQ, double muS){
   return Chi2STaylor(T,muB,muS,muQ);
}

double P2BQ(double T, double muB, double muQ, double muS){
   return Chi11BQTaylor(T,muB,muS,muQ);
}

double P2BS(double T, double muB, double muQ, double muS){
   return Chi11BSTaylor(T,muB,muS,muQ);
}

double P2QS(double T, double muB, double muQ, double muS){
   return Chi11QSTaylor(T,muB,muS,muQ);
}

double P2TB(double T, double muB, double muQ, double muS){
   return DBarDensDTTaylor(T,muB,muS,muQ);
}

double P2TQ(double T, double muB, double muQ, double muS){
   return DChDensDTTaylor(T,muB,muS,muQ);
}

double P2TS(double T, double muB, double muQ, double muS){
   return DStrDensDTTaylor(T,muB,muS,muQ);
}

double P2T2(double T, double muB, double muQ, double muS){
   return DEntrDTTaylor(T,muB,muS,muQ);
}
#undef NRANSI
