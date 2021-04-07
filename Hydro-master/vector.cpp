#ifndef _VECTOR_CPP_
#define _VECTOR_CPP_

#include <math.h>
#include "vector.h"

//not sure if need this???
template <class T, int D>
VectorConsts<T,D>::VectorConsts() {

	for(int i=0; i<D; i++) {
		uni[i].x[i]=(T)1;
		Uni.x[i]=(T)1;
	}
	
}

template <int D>
void unitaryVector(Vector<double,D> & V) 
{ 
	for(int i=0; i<D; i++) V.x[i]=1.0;
}

template <int D>
void unitaryVector(Vector<int,D> & V) 
{ 
	for(int i=0; i<D; i++) V.x[i]=1;
}


template <class T, int D>
Vector<T,D>::Vector() {
                      
	for(int i=0; i<D; i++) x[i]=0;
	
}


template <class T, int D>
template <class U>
Vector<T,D>& Vector<T,D>::operator=(Vector<U,D> a) {
             
	for(int i=0; i<D; i++) x[i]=(T)a.x[i];
	
	return *this;
	
}





template <class T, int D>
Vector<T,D>& Vector<T,D>::operator=(double a) {
             
	for(int i=0; i<D; i++) x[i]=(T)a;
	
	return *this;
	
}


template <class T, int D>
Vector<T,D>& Vector<T,D>::operator+=(Vector<T,D> a) {
             
	for(int i=0; i<D; i++) x[i]+=a.x[i];
	
	return *this;
	
}


template <class T, int D>
Vector<T,D>& Vector<T,D>::operator-=(Vector<T,D> a) {
             
	for(int i=0; i<D; i++) x[i]-=a.x[i];
	
	return *this;
	
}


template <class T, int D>
Vector<T,D>& Vector<T,D>::operator*=(T l) {
             
	int i;
	
	for(i=0; i<D; i++) x[i]*=l;
	
	return *this;
	
}


template <class T, int D>
Vector<T,D> operator+ (Vector<T,D> a, Vector<T,D> b) {
            
	Vector<T,D> t;
	
	t=0;
	
	return (t+=a)+=b;
	
}


template <class T, int D>
Vector<T,D> operator- (Vector<T,D> a) {
            
	Vector<T,D> t;
	
	t=0;
	
	return t-=a;
	
}


template <class T, int D>
Vector<T,D> operator- (Vector<T,D> a, Vector<T,D> b) {
            
	Vector<T,D> t;
	
	return t=a+(-b);
	
}


template <class T, int D>
Vector<T,D> operator* (T l, Vector<T,D> a) {
            
	Vector<T,D> t;
	
	return (t=a)*=l;
	
}


template <class T, int D>
Vector<T,D> Direction (Vector<T,D> a) {
            
	double t=Norm(a);
	
	if(t==0) {Vector<T,D> u; return u;}
	
	return 1/Norm(a)*a;
	
}


template <class T, int D>
double Norm(Vector<T,D> a) {
       
	int i;
	double t=0;
	
	for(i=0; i<D; i++) t+=a.x[i]*a.x[i];
	
	return sqrt(t);
	
}


template <class T, int D>
double Norm2(Vector<T,D> a) {
       
	double t=0;
	
	for(int i=0; i<D; i++) t+=a.x[i]*a.x[i];
	
	return t;
	
}


template <class T, int D>
ostream& operator<<(ostream& os, Vector<T,D> a) {
         
	for(int i=0; i<D; i++) os<<a.x[i]<<" ";
	
	return os;
	
}


template <class T, int D> 
Vector<T,D> zeroVector(void) {
            
	Vector<T,D> t;
	
	for(int i=0; i<D; i++) t.x[i]=0;
	
	return t;
	
}


template <class T, int D>
Vector<T,D> unitaryVector(int i) {
            
	if( i>=D ) return zeroVector<T,D>();
	
	Vector<T,D> t=zeroVector<T,D>();
	
	t.x[i]=1;
	
	return t;
	
}


template <class T, int D>
Vector<T,D> unitaryVector() {
            
	Vector<T,D> t=zeroVector<T,D>();
	
	for(int i=0; i<D; i++) t.x[i]=1;
	return t;
	
}


template <class T, int D> 
double inner (Vector<T,D> a, Vector<T,D> b) {
       
	double t=0;
	
	for(int i=0; i<D; i++) t+=a.x[i]*b.x[i];
	
	return t;
	
}


#endif
