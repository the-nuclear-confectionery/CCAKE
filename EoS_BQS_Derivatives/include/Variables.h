//  Variables.h

#ifndef Variables_h
#define Variables_h

/* For parameters inpu8t. */
double C0, A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, B0, B1, B2, B3, B4, B5, B6, B7, B8, B9;

double **parMatrix, *par;

int i,j,k,l;

double *xSN, *fSN;
int check;

double Tval, muBval, muSval, muQval;
double Rat;

double PressVal, EntrVal, BarDensVal, StrDensVal, ChDensVal, EnerDensVal, Chi3BVal, IsEntrVal, SpSoundVal;

double *CHI000PAR,*CHI200PAR,*CHI020PAR,*CHI002PAR,*CHI110PAR,*CHI101PAR,*CHI011PAR,*CHI400PAR,*CHI040PAR,*CHI004PAR,*CHI310PAR,*CHI301PAR,
        *CHI031PAR,*CHI130PAR,*CHI103PAR,*CHI013PAR,*CHI220PAR,*CHI202PAR,*CHI022PAR,*CHI211PAR,*CHI121PAR,*CHI112PAR;
        
double C2B2, C2Q2, C2S2, C2BQ, C2BS, C2QS, C2TB, C2TQ, C2TS, C2T2;
double C1B, C1Q, C1S, C1T;
double  D2PB2, D2PQ2, D2PS2, D2PBQ, D2PBS, D2PQS, D2PTB, D2PTQ, D2PTS, D2PT2;

#endif /* Variables_h */
