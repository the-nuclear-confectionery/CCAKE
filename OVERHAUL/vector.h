#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include "mathdef.h"
#include "vector.h"

using namespace std;

template <class T, int D>
class Vector {
public:
    T x[D];
    Vector<T,D>();
    template <class U> Vector<T,D>& operator=(Vector<U,D>);
    Vector<T,D>& operator=(double);
    Vector<T,D>& operator+=(Vector<T,D>);
    Vector<T,D>& operator-=(Vector<T,D>);
    Vector<T,D>& operator*=(T);
};

template <class T, int D>
class VectorConsts {
public:
    Vector<T,D> uni[D];
    Vector<T,D> Uni;
    VectorConsts<T,D>();
};

typedef Vector<double,1> vec1;
typedef Vector<double,2> vec2;
typedef Vector<double,3> vec3;
typedef Vector<int,1> tri1;
typedef Vector<int,2> tri2;
typedef Vector<int,3> tri3;

template <class T, int D> Vector<T,D> operator+ (Vector<T,D>, Vector<T,D>);
template <class T, int D> Vector<T,D> operator- (Vector<T,D>);
template <class T, int D> Vector<T,D> operator- (Vector<T,D>, Vector<T,D>);
template <class T, int D> Vector<T,D> operator* (T, Vector<T,D>);
template <class T, int D> double inner (Vector<T,D>, Vector<T,D>);
template <class T, int D> double Norm(Vector<T,D>);
template <class T, int D> double Norm2(Vector<T,D>);
template <class T, int D> Vector<T,D> Direction(Vector<T,D>);
template <class T, int D> ostream& operator<<(ostream&, Vector<T,D>);
template <int D> void unitaryVector(Vector<double,D> & V);
template <int D> void unitaryVector(Vector<int,D> & V);


#include "vector.cpp"

#endif
