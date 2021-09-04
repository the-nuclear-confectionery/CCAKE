OBJS = main.o vector.o hydrosim.o newspline.o random.o eos.o matrix.o enteric.o
#CC = g++
CC = h5c++
DEBUG = -g
STANDARD= -std=c++11
WARNINGS= -w -Wall -Wextra -pedantic -Wno-conversion -Wcast-align -Wunused -Wshadow -Wpointer-arith
OPTIMIZATION= -Og
#OPTIMIZATION= -Ofast -funroll-loops -fvectorize -finline-functions -s -DNDEBUG -march=native -mtune=native
# change to -O3 or -Ofast for production
EIGENPATH=$(realpath ../eigen)
SPLINTERPATH=$(realpath ../splinter/include)
SPLINTERBUILDPATH=$(realpath ../splinter/build)
GSLPATH=/projects/jnorhos/local/include
GSLLIBPATH=/projects/jnorhos/local/lib
EOSLIBPATH=$(realpath ../EoS_BQS_Derivatives/lib)
CFLAGS = $(WARNINGS) $(OPTIMIZATION) -c -I$(EOSLIBPATH) -I$(EIGENPATH) -I$(SPLINTERPATH) -I$(GSLPATH) $(DEBUG) $(STANDARD)
LDFLAGS = -L$(EOSLIBPATH) -lbsqeos -lm -lgsl -lgslcblas -lblas -lstdc++ -L$(GSLLIBPATH) -L$(SPLINTERBUILDPATH) -lsplinter-3-0 $(DEBUG) $(STANDARD)

export HDF5_LDFLAGS+=$(LDFLAGS)
export HDF5_LIBS+=-L$(SPLINTERBUILDPATH) -lsplinter-3-0
#export HDF5_CXXFLAGS=-I$(SPLINTERPATH)


vusphydro : $(OBJS)
	#$(info [${LD_RUN_PATH}])
	#$(info [${LD_LIBRARY_PATH}])
	#$(info [${HDF5_LDFLAGS}])
	#$(info [${HDF5_LIBS}])
	#$(info [${HDF5_CXXFLAGS}])
	$(CC) -g $(OBJS) -o vusphydro $(LDFLAGS)

main.o: vector.h mathdef.h particle.h eos.h enteric.h LinkList.h random.h matrix.h hydrosim.h eostables.h read_in_hdf/read_in_hdf.h main.cpp
	$(CC) $(CFLAGS) main.cpp

vector.o: vector.h mathdef.h vector.cpp
	$(CC) $(CFLAGS) vector.cpp

matrix.o: matrix.cpp matrix.h vector.h mathdef.h
	$(CC) $(CFLAGS) matrix.cpp

hydrosim.o: vector.h mathdef.h particle.h hydrosim.h LinkList.h rungekutta4.h eos.h output.h bbmg.h matrix.h hydrosim.cpp
	$(CC) $(CFLAGS) hydrosim.cpp

newspline.o: newspline.cpp newspline.hpp
	$(CC) $(CFLAGS) newspline.cpp

eos.o: tables.h eos.h eostables.h eos.cpp $(SPLINTERPATH)/datatable.h $(SPLINTERPATH)/bspline.h $(SPLINTERPATH)/bsplinebuilder.h read_in_hdf/read_in_hdf.h
	$(CC) $(CFLAGS) eos.cpp

random.o: random.cpp random.h
	$(CC) $(CFLAGS) random.cpp

enteric.o:  vector.h particle.h eos.h tables.h enteric.h LinkList.h matrix.h eostables.h enteric.cpp
	$(CC) $(CFLAGS) enteric.cpp



# deletes all .o files in case of error in makefile.  This should not be something needed regularly, only if you are updating the code and run into an error.
clean:
	-rm *.o vusphydro

# puts all the files in a vusphydro.tar file
#tar:
#    tar cfv vusphydro.tar vector.h main.cpp vector.cpp
