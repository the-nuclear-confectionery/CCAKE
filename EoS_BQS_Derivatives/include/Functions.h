#ifndef Functions_h
#define Functions_h

double **parMatrix;
double *par;

void set_lowT_parameters(double *par1, double *par2);
void set_lowT_Mod_parameters(double *par1, double *par2);


// ------------- THE FUNCTIONAL FORMS --------------- //
double coeff(double *par,double x); 
double coeffprime(double *par,double x);

// ------------- THE COEFFICIENTS ------------------ //
double CHI000(double T);
double CHI200(double T);
double CHI020(double T);
double CHI002(double T);
double CHI110(double T);
double CHI101(double T);
double CHI011(double T);
double CHI400(double T);
double CHI040(double T);
double CHI004(double T);
double CHI310(double T);
double CHI301(double T);
double CHI031(double T);
double CHI130(double T);
double CHI103(double T);
double CHI013(double T);
double CHI220(double T);
double CHI202(double T);
double CHI022(double T);
double CHI211(double T);
double CHI121(double T);
double CHI112(double T);

// ---- Derivatives wrt T --------- //
double DCHI000DT(double T);
double DCHI200DT(double T);
double DCHI020DT(double T);
double DCHI002DT(double T);
double DCHI110DT(double T);
double DCHI101DT(double T);
double DCHI011DT(double T);
double DCHI400DT(double T);
double DCHI040DT(double T);
double DCHI004DT(double T);
double DCHI310DT(double T);
double DCHI301DT(double T);
double DCHI031DT(double T);
double DCHI130DT(double T);
double DCHI103DT(double T);
double DCHI013DT(double T);
double DCHI220DT(double T);
double DCHI202DT(double T);
double DCHI022DT(double T);
double DCHI211DT(double T);
double DCHI121DT(double T);
double DCHI112DT(double T);

// ---- Double Derivatives wrt T --------- //
double D2CHI000DT2(double T);
double D2CHI200DT2(double T);
double D2CHI020DT2(double T);
double D2CHI002DT2(double T);
double D2CHI110DT2(double T);
double D2CHI101DT2(double T);
double D2CHI011DT2(double T);
double D2CHI400DT2(double T);
double D2CHI040DT2(double T);
double D2CHI004DT2(double T);
double D2CHI310DT2(double T);
double D2CHI301DT2(double T);
double D2CHI031DT2(double T);
double D2CHI130DT2(double T);
double D2CHI103DT2(double T);
double D2CHI013DT2(double T);
double D2CHI220DT2(double T);
double D2CHI202DT2(double T);
double D2CHI022DT2(double T);
double D2CHI211DT2(double T);
double D2CHI121DT2(double T);
double D2CHI112DT2(double T);

double lowT_CHI000(double T);
double lowT_CHI200(double T);
double lowT_CHI020(double T);
double lowT_CHI002(double T);
double lowT_CHI110(double T);
double lowT_CHI101(double T);
double lowT_CHI011(double T);
double lowT_CHI400(double T);
double lowT_CHI040(double T);
double lowT_CHI004(double T);
double lowT_CHI310(double T);
double lowT_CHI301(double T);
double lowT_CHI031(double T);
double lowT_CHI130(double T);
double lowT_CHI103(double T);
double lowT_CHI013(double T);
double lowT_CHI220(double T);
double lowT_CHI202(double T);
double lowT_CHI022(double T);
double lowT_CHI211(double T);
double lowT_CHI121(double T);
double lowT_CHI112(double T);

double lowT_DCHI000DT(double T);
double lowT_DCHI200DT(double T);
double lowT_DCHI020DT(double T);
double lowT_DCHI002DT(double T);
double lowT_DCHI110DT(double T);
double lowT_DCHI101DT(double T);
double lowT_DCHI011DT(double T);
double lowT_DCHI400DT(double T);
double lowT_DCHI040DT(double T);
double lowT_DCHI004DT(double T);
double lowT_DCHI310DT(double T);
double lowT_DCHI301DT(double T);
double lowT_DCHI031DT(double T);
double lowT_DCHI130DT(double T);
double lowT_DCHI103DT(double T);
double lowT_DCHI013DT(double T);
double lowT_DCHI220DT(double T);
double lowT_DCHI202DT(double T);
double lowT_DCHI022DT(double T);
double lowT_DCHI211DT(double T);
double lowT_DCHI121DT(double T);
double lowT_DCHI112DT(double T);

double lowT_D2CHI000DT2(double T);
double lowT_D2CHI200DT2(double T);
double lowT_D2CHI020DT2(double T);
double lowT_D2CHI002DT2(double T);
double lowT_D2CHI110DT2(double T);
double lowT_D2CHI101DT2(double T);
double lowT_D2CHI011DT2(double T);
double lowT_D2CHI400DT2(double T);
double lowT_D2CHI040DT2(double T);
double lowT_D2CHI004DT2(double T);
double lowT_D2CHI310DT2(double T);
double lowT_D2CHI301DT2(double T);
double lowT_D2CHI031DT2(double T);
double lowT_D2CHI130DT2(double T);
double lowT_D2CHI103DT2(double T);
double lowT_D2CHI013DT2(double T);
double lowT_D2CHI220DT2(double T);
double lowT_D2CHI202DT2(double T);
double lowT_D2CHI022DT2(double T);
double lowT_D2CHI211DT2(double T);
double lowT_D2CHI121DT2(double T);
double lowT_D2CHI112DT2(double T);


// ------------ THERMODYNAMICS --------------------- //
double PressTaylor(double T, double muB, double muQ, double muS);
double EntrTaylor(double T, double muB, double muQ, double muS);
double BarDensTaylor(double T, double muB, double muQ, double muS);
double StrDensTaylor(double T, double muB, double muQ, double muS);
double ChDensTaylor(double T, double muB, double muQ, double muS);

double Chi2BTaylor(double T, double muB, double muQ, double muS);
double Chi2QTaylor(double T, double muB, double muQ, double muS);
double Chi2STaylor(double T, double muB, double muQ, double muS);
double Chi11BQTaylor(double T, double muB, double muQ, double muS);
double Chi11BSTaylor(double T, double muB, double muQ, double muS);
double Chi11QSTaylor(double T, double muB, double muQ, double muS);

double DBarDensDTTaylor(double T, double muB, double muQ, double muS);
double DChDensDTTaylor(double T, double muB, double muQ, double muS);
double DStrDensDTTaylor(double T, double muB, double muQ, double muS);

double DEntrDTTaylor(double T, double muB, double muQ, double muS);

double SpSound(double T, double muB, double muQ, double muS);
double SpSoundReadable(double T, double muB, double muQ, double muS);

double P2B2(double T, double muB, double muQ, double muS);
double P2Q2(double T, double muB, double muQ, double muS);
double P2S2(double T, double muB, double muQ, double muS);
double P2BQ(double T, double muB, double muQ, double muS);
double P2BS(double T, double muB, double muQ, double muS);
double P2QS(double T, double muB, double muQ, double muS);
double P2TB(double T, double muB, double muQ, double muS);
double P2TQ(double T, double muB, double muQ, double muS);
double P2TS(double T, double muB, double muQ, double muS);
double P2T2(double T, double muB, double muQ, double muS);

#endif /* Functions_h */
