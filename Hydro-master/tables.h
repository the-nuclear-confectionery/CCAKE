#ifndef TABLES_H_
#define TABLES_H_

#include <string>
#include "vector.h"
#include "matrix.h"

typedef struct _Table1 {

    double mub,mus;                   // entropy density
    double T,s;                   // temperature - used for entropy table
} Table1;

typedef struct _Table2 {

    double p;                   // pressure
    double e;                   // energy density
    double s;              // derivative of entalpy in entropy
    double nb;

} Table2;

struct _inputIC {

    int on,start,end; // on=1 if input of start and end events are given at execution
    string man,rnum;
    double dt;
};

extern char ifolder [];
extern std::string ofolder;
extern double  freezeoutT,zconst;
extern  Matrix <double,2,2> Imat;
extern int hcor;

#endif
