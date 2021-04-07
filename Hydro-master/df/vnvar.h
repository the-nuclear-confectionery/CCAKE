#ifndef _VNVAR_H_
#define _VNVAR_H_

#include "vector.h"

class VNVAR
{
private:
	int ptmax;

public:
	Vector<double,2> I,Q;
	double psi,intert;
	VNVAR();
	void restart();

};

VNVAR::VNVAR()
{
	psi=0;
}





void VNVAR::restart()
{
	I=0;
	Q=0;
	intert=0;
	psi=0;
	
}

#endif
