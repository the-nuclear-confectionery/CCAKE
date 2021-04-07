#ifndef BESSEL_H_ 
#define BESSEL_H_


using namespace std;


class Bessel {
public:
	double I0(double x);
	double I1(double x);
	double K0(double x);
	double K1(double x);
	double Kn(int n, double x);
};

#endif


