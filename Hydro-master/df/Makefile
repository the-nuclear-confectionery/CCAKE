OBJS = main.o bessel.o
CC = c++
DEBUG = -g 
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

sphrio : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o fo

main.o: SPH.h vector.h int.h tables.h vnvar.h spectra.h bessel.h main.cpp 
	$(CC) $(CFLAGS) main.cpp

bessel.o: bessel.h bessel.cpp 
	$(CC) $(CFLAGS) bessel.cpp
	

# deletes all .o files in case of error in makefile.  This should not be something needed regularly, only if you are updating the code and run into an error.	
#clean:
#	\rm *.o *~ sphrio

# puts all the files in a sphrio.tar file
#tar:
#    tar cfv sphrio.tar vector.h main.cpp vector.cpp 
