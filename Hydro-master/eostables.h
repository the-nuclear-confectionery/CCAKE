#ifndef EOSTABLES_H_ 
#define EOSTABLES_H_



typedef struct _eostableshigh {
     
	double T,s,p,e,dtds;               // temperature - used for entropy table	
	int low;
} eostableshigh;

typedef struct _eostableslow {
     
	double T,s,p,e,dtds;	
    
} eostableslow;


extern eostableshigh *ETH;
extern eostableslow *ETL;
extern int nETH,nETL;

#endif
