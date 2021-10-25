#ifndef MATRIX_CPP
#define MATRIX_CPP

#include <math.h>
#include "matrix.h"
#include "vector.h"
#include <stdio.h>

//not sure if need this???

// D1 is the number for the rows, D2 is the number for the columns
template <class T, int D1, int D2>
Matrix<T,D1,D2>::Matrix() {

    for(int i=0; i<D1; i++) {
        for(int j=0; j<D2; j++) {
            x[i][j]=0;
        }
    }

}


template <class T, int D1, int D2>
template <class U>
Matrix<T,D1,D2>& Matrix<T,D1,D2>::operator=(Matrix<U,D1,D2> a) {

    for(int i=0; i<D1; i++) {
        for(int j=0; j<D2; j++) {
            x[i][j]=(T)a.x[i][j];
        }
    }

    return *this;

}





template <class T, int D1, int D2>
Matrix<T,D1,D2>& Matrix<T,D1,D2>::operator=(double a) {

    for(int i=0; i<D1; i++) {
        for(int j=0; j<D2; j++) {
            x[i][j]=(T)a;
        }
    }

    return *this;

}


template <class T, int D1, int D2>
Matrix<T,D1,D2>& Matrix<T,D1,D2>::operator+=(Matrix<T,D1,D2> a) {

    for(int i=0; i<D1; i++) {
        for(int j=0; j<D2; j++) {
            x[i][j]+=a.x[i][j];
        }
    }

    return *this;

}


template <class T, int D1, int D2>
Matrix<T,D1,D2>& Matrix<T,D1,D2>::operator-=(Matrix<T,D1,D2> a) {

    for(int i=0; i<D1; i++) {
        for(int j=0; j<D2; j++) {
            x[i][j]-=a.x[i][j];
        }
    }

    return *this;

}


template <class T, int D1, int D2>
Matrix<T,D1,D2>& Matrix<T,D1,D2>::operator*=(T l) {



    for(int i=0; i<D1; i++) {
        for(int j=0; j<D2; j++) {
            x[i][j]*=l;
        }
    }

    return *this;

}

//only works for matrices with the same dimensions
template <class T, int D1, int D2>
Matrix<T,D1,D2>& Matrix<T,D1,D2>::operator*=(Matrix<T,D2,D1> b)
{


    for(int i=0; i<D1; i++) {
        for(int j=0; j<D1; j++) {
            double sub=0;
            for(int k=0; k<D2; k++) {
                sub+=x[i][k]*b.x[k][j];
            }
            x[i][j]=sub;

        }
    }

    return *this;

}


template <class T, int D1, int D2>
Matrix<T,D1,D2> operator+ (Matrix<T,D1,D2> a, Matrix<T,D1,D2> b) {

    Matrix<T,D1, D2> t;

    t=0;

    return (t+=a)+=b;

}


template <class T, int D1, int D2>
Matrix<T,D1,D2> operator- (Matrix<T,D1,D2> a) {

    Matrix<T,D1,D2> t;

    t=0;

    return t-=a;

}


template <class T, int D1, int D2>
Matrix<T,D1,D2> operator- (Matrix<T,D1,D2> a, Matrix<T,D1,D2> b) {

    Matrix<T,D1,D2> t;

    return t=a+(-b);

}


template <class T, int D1, int D2>
Matrix<T,D1,D2> operator* (T l, Matrix<T,D1,D2> a) {

    Matrix<T,D1,D2> t;

    return (t=a)*=l;

}

template <class T, int D1, int D2>
Matrix<T,D1,D2>& Matrix<T,D1,D2>::identity() {

    if (D1!=D2) cout << "Error: not true indentity matrix!" << endl;


    for(int i=0; i<D1; i++) {
        for(int j=0; j<D1; j++) {
            if (i==j) x[i][j]=1;
            else x[i][j]=0;
        }
    }

    return *this;

}



template <class T, int D1, int D2>
ostream& operator<<(ostream& os, Matrix<T,  D1, D2> a) {

    for(int j=0; j<D2; j++) {
        for(int i=0; i<D1; i++) {
            os << a.x[i][j] <<" ";
        }
    }

    return os;

}

template <class T, int D1, int D2>
Matrix<T,D2,D1> transpose(Matrix<T,D1,D2> a)
{
    Matrix<T,D1,D2> t;


    for(int i=0; i<D2; i++) {
        for(int j=0; j<D1; j++) {
            t.x[i][j]=a.x[j][i];
        }
    }

    return t;
}

template <class T, int D1, int D2, int Da2, int Db1>
Matrix<T,D1,D2> operator* (Matrix<T, D1,Da2> a, Matrix< T,Db1, D2> b) {

    if (Da2!=Db1) {
        cout << "Error: Attempt to multiple two matrices of improper sizes!" << endl;
        getchar();
    }


    Matrix<T,D1, D2> t;

    for(int i=0; i<D2; i++) {
        for(int j=0; j<D1; j++) {
            t.x[j][i]=0;
            for(int k=0; k<Da2; k++) {
                t.x[j][i]+=a.x[j][k]*b.x[k][i];
            }


        }
    }

    return t;

}

template <class T, int D1, int D2>
Vector< T,D1> operator* (Matrix<T, D1,D2> a, Vector< T,D2> b) {

    Vector< T,D1> t;


    for(int j=0; j<D1; j++) {
        t.x[j]=0;
        for(int k=0; k<D2; k++) {
            t.x[j]+=a.x[j][k]*b.x[k];
        }
    }

    return t;

}


template <class T, int D1, int D2>
Matrix<T, D1, D2> operator* (Vector< T,D1> a, Vector<T,D2> b) {




    Matrix<T,  D1,D2> t;

    for(int i=0; i<D2; i++) {
        for(int j=0; j<D1; j++) {
            t.x[i][j]=a.x[i]*b.x[j];


        }
    }

    return t;

}

template <class T, int D1, int D2>
Vector<T,D2> column(int l, Matrix<T, D1, D2> a)
{
    Vector<T,D2> v;

    for(int i=0; i<D1; i++) {
        v.x[i]=a.x[i][l];

    }

    return v;


}

template <class T, int D1, int D2>
Vector<T,D1> row(int l, Matrix<T, D1,D2> a)
{
    Vector<T,D1> v;

    for(int i=0; i<D2; i++) {
        v.x[i]=(T)a.x[l][i];

    }

    return v;


}

template <class T, int D1, int D2>
Vector<T,(D2-1)> rowm1(int l, Matrix<T, D1,D2> a)
{
    Vector<T,(D2-1)> v;

    for(int i=0; i<(D2-1); i++) {
        v.x[i]=(T)a.x[l][i];

    }

    return v;


}



template <class T, int D1, int D2>
Vector<T,(D2-1)> rowp1(int l, Matrix<T, D1,D2> a)
{
    Vector<T,(D2-1)> v;

    for(int i=1; i<D2; i++) {
        v.x[i-1]=(T)a.x[l][i];

    }




    return v;


}

template <class T, int D1, int D2>
Vector<T,(D1-1)> colp1(int l, Matrix<T, D1,D2> a)
{
    Vector<T,(D1-1)> v;

    for(int i=1; i<D1; i++) {
        v.x[i-1]=(T)a.x[i][l];

    }

    return v;


}

template <class T, int D1, int D2>
void  mini(Matrix<T, D1-1, D2-1> &b, Matrix<T, D1, D2> a)
{


    for(int j=1; j<=2; j++) {
        for(int i=1; i<=2; i++) {
            b.x[i-1][j-1]=(T)a.x[i][j];

        }
    }



}

template <class T, int D1, int D2>
Vector< T,D2> operator* ( Vector< T,D1> a, Matrix<T, D1,D2> b) {

    Vector< T,D2> t;


    for(int j=0; j<D2; j++) {
        t.x[j]=0;
        for(int k=0; k<D1; k++) {
            t.x[j]+=a.x[k] *b.x[k][j];
        }
    }

    return t;

}

template <class T,  int D2>
Matrix<T, 2, D2> stack(Vector<T, D2> a)
{
    Matrix<T,2,D2> t;


    for(int i=0; i<2; i++) {
        for(int j=0; j<D2; j++) {
            t.x[i][j]=a.x[j];
        }
    }

    return t;
}

template <class T, int D1>
Matrix<T,D1,D1> diagonal(Vector<T, D1> a)
{
    Matrix<T,D1,D1> t;

    for(int i=0; i<D1; i++) {
        for(int j=0; j<D1; j++) {
            if (i==j) t.x[i][j]=a.x[i];
            else t.x[i][j]=0;
        }
    }

    return t;

}

template <class T, int D1> double deter(Matrix<T, D1, D1> a)
{
    double d;
    d=a.x[0][0]*a.x[1][1]-a.x[0][1]*a.x[1][0];

    return d;
}

template <class T, int D1, int D2>
double con(Matrix<T, D1,D2> a, Matrix< T,D1, D2> b) {




    double t=0;

    for(int i=0; i<D2; i++) {
        for(int j=0; j<D1; j++) {
            t+=a.x[j][i]*b.x[j][i];
        }
    }

    return t;

}

template <class T, int D1, int D2>
double con2(Matrix<T, D1,D2> a, Matrix< T,D1, D2> b) {




    double t=0.;

    for(int i=0; i<D2; i++) {
        for(int j=0; j<D1; j++) {
            t+=a.x[i][j]*b.x[i][j];
        }
    }

    return t;

}

template <class T, int D1, int D2>
void tmini( Matrix<T, D1, D2> &b, Matrix<T, D1-1, D2-1>a)
{



    for(int j=0; j<=(D1-1); j++) {
        for(int i=0; i<=(D1-1); i++) {
            b.x[i+1][j+1]=(T)a.x[i][j];

        }
    }



}


#endif
