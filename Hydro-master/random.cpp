#ifndef _RANDOM_CPP_
#define _RANDOM_CPP_

#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>
#include "random.h"

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <vector>

using namespace std;



double randall(double min, double max, int & r) {

    srand( time( NULL ) +r);


    double out=min+ rand()* (max-min) / double(RAND_MAX+1.);

    r=r+rand();

    return out;

} // random for any min/max, doubles only

double rand0(double max,int & r) {

    srand( time( NULL )+r );

    double out=rand()*( max) / double(RAND_MAX+1.);

    r=r+rand();
    return out;

} // random for min=0, doubles only

double rand01(int & r) {

    srand( time( NULL )+r );

    double out=rand() / double(RAND_MAX+1.);

    r=r+rand();
    return out;

} // random for min=0, max=1, doubles only


int randint(int max,int & r ) {

    srand( time( NULL )+r );


    int out= (rand() % max);


    r=r+rand();

    return out;
}



#endif
