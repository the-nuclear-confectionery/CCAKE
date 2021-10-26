OBJS = main.o vector.o matrix.o BSQHydro.o runge_kutta.o particle.o eos.o system_state.o equations_of_motion.o eos_delaunay.o
#CC = g++
CC = h5c++
DEBUG = -g
STANDARD= -std=c++11
WARNINGS= -w -Wall -Wextra -pedantic -Wno-conversion -Wcast-align -Wunused -Wshadow -Wpointer-arith
OPTIMIZATION= -Og
#OPTIMIZATION= -Ofast -funroll-loops -fvectorize -finline-functions -s -DNDEBUG -march=native -mtune=native
# change to -O3 or -Ofast for production
EIGENPATH=$(realpath ../eigen)
GSLPATH=/projects/jnorhos/local/include
GSLLIBPATH=/projects/jnorhos/local/lib
EOSLIBPATH=$(realpath ../EoS_BQS_Derivatives/lib)
CFLAGS = $(WARNINGS) $(OPTIMIZATION) -c -I$(EOSLIBPATH) -I$(EIGENPATH) \
        -I$(GSLPATH) $(DEBUG) $(STANDARD)
LDFLAGS = -L$(EOSLIBPATH) -lbsqeos -lm -lgsl -lgslcblas -lblas -lstdc++ \
          -L$(GSLLIBPATH) $(DEBUG) $(STANDARD)

export HDF5_LDFLAGS+=$(LDFLAGS)

vusphydro : $(OBJS)
	$(CC) -g $(OBJS) -o vusphydro $(LDFLAGS)

main.o: vector.h mathdef.h particle.h eos.h input_output.h linklist.h matrix.h BSQHydro.h read_in_hdf/read_in_hdf.h main.cpp
	$(CC) $(CFLAGS) main.cpp

vector.o: vector.h mathdef.h vector.cpp
	$(CC) $(CFLAGS) vector.cpp

matrix.o: matrix.cpp matrix.h vector.h mathdef.h
	$(CC) $(CFLAGS) matrix.cpp

runge_kutta.o: runge_kutta.cpp runge_kutta.h equations_of_motion.h system_state.h
	$(CC) $(CFLAGS) runge_kutta.cpp

system_state.o: system_state.cpp system_state.h
	$(CC) $(CFLAGS) system_state.cpp


BSQHydro.o: vector.h mathdef.h particle.h BSQHydro.h linklist.h runge_kutta.h eos.h matrix.h BSQHydro.cpp
	$(CC) $(CFLAGS) BSQHydro.cpp

particle.o: particle.cpp particle.h eos.h
	$(CC) $(CFLAGS) particle.cpp

eos.o: eos.h eos.cpp read_in_hdf/read_in_hdf.h eos_delaunay/eos_delaunay.h eos_delaunay/eos_delaunay.cpp
	$(CC) $(CFLAGS) eos.cpp

equations_of_motion.o:  vector.h particle.h eos.h equations_of_motion.h
	$(CC) $(CFLAGS) equations_of_motion.cpp

eos_delaunay.o: eos_delaunay/eos_delaunay.cpp eos_delaunay/eos_delaunay.h
	$(CC) $(CFLAGS) $(INCFLAGS) -fpic -c eos_delaunay/eos_delaunay.cpp 
	
	

clean:
	-rm *.o vusphydro
