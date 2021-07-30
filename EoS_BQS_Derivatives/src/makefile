EoS_BQS:
	#gcc -o EoS_BQS *.c -lm
	h5cc -o EoS_BQS *.c -I/projects/jnorhos/local/gsl/include -L/projects/jnorhos/local/gsl/lib -lm -lgsl -lgslcblas -lblas

install:
	cp EoS_BQS ../
	
clean:
	rm EoS_BQS
