EoS_BQS:
	#gcc -o EoS_BQS *.c -lm
	h5cc -o EoS_BQS *.c -I$(GSLPATH) -L$(GSLLIBPATH) -lm -lgsl -lgslcblas -lblas

install:
	cp EoS_BQS ../
	
clean:
	rm EoS_BQS
