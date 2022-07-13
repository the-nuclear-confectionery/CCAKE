FLAGS_NEEDED_FOR_GCC_10_3_0 = -fcommon

EoS_BQS:
	#gcc -o EoS_BQS *.c -lm
	h5cc $(FLAGS_NEEDED_FOR_GCC_10_3_0) -o EoS_BQS *.c -I/projects/jnorhos/local/include -L/projects/jnorhos/local/lib -lm -lgsl -lgslcblas -lblas

install:
	cp EoS_BQS ../
	
clean:
	rm -f EoS_BQS *.o
