EoS_BQS:
	GSLPATH=/projects/jnorhos/local/gsl/include
	GSLLIBPATH=/projects/jnorhos/local/gsl/lib
	#gcc -o EoS_BQS *.c -lm
	h5cc -o EoS_BQS *.c -I$(GSLPATH) -L$(GSLLIBPATH) -lm -lgsl -lgslcblas -lblas

install:
	cp EoS_BQS ../
	
clean:
	rm EoS_BQS
