OBJDIR=obj
OBJS = $(OBJDIR)/transport_coefficients.o $(OBJDIR)/main.o $(OBJDIR)/vector.o \
			 $(OBJDIR)/matrix.o $(OBJDIR)/BSQHydro.o $(OBJDIR)/particle.o $(OBJDIR)/eos.o \
       $(OBJDIR)/system_state.o $(OBJDIR)/linklist.o $(OBJDIR)/eos_delaunay.o \
			 $(OBJDIR)/kernel.o $(OBJDIR)/sph_workstation.o $(OBJDIR)/input_output.o \
			 $(OBJDIR)/rootfinder.o $(OBJDIR)/eos_table.o $(OBJDIR)/eos_extension.o
			 
#CC = g++
CC = h5c++

#DEBUG = -g -D_GLIBCXX_CONCEPT_CHECKS -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
DEBUG=-g
STANDARD= -std=c++11
WARNINGS= -w -Wall -Wextra -pedantic -Wno-conversion -Wcast-align -Wunused \
          -Wshadow -Wpointer-arith
OPTIMIZATION= -Og
#OPTIMIZATION= -Ofast -s -march=native -mtune=native
# change to -O3 or -Ofast for production

SRC=$(realpath src)
INCLUDE=$(realpath include)

#EIGENPATH=$(realpath ../eigen)
GSLPATH=/projects/jnorhos/local/include
GSLLIBPATH=/projects/jnorhos/local/lib
EOSLIBPATH=$(realpath ../EoS_BQS_Derivatives/lib)
CFLAGS = $(WARNINGS) $(OPTIMIZATION) -c -I$(EOSLIBPATH) -I$(EIGENPATH) \
        -I$(GSLPATH) $(DEBUG) $(STANDARD) -fopenmp
LDFLAGS = -L$(EOSLIBPATH) -lbsqeos -lm -lgsl -lgslcblas -lblas -lstdc++ \
          -L$(GSLLIBPATH) $(DEBUG) $(STANDARD) -fopenmp

#export OMP_NUM_THREADS=`nproc --all`
export OMP_NUM_THREADS=1
export HDF5_LDFLAGS+=$(LDFLAGS)

.PHONY:	all clean mkobjdir persephone

all:	mkobjdir persephone
	
#target:	mkobjdir persephone

persephone: mkobjdir $(OBJS)
	(cd $(EOSLIBPATH); make clean; make)
	$(CC) -g $(OBJS) -o persephone $(LDFLAGS)

mkobjdir:
	-@mkdir -p $(OBJDIR)

$(OBJDIR)/main.o: $(INCLUDE)/vector.h $(INCLUDE)/mathdef.h $(INCLUDE)/particle.h \
        $(INCLUDE)/eos.h $(INCLUDE)/input_output.h $(INCLUDE)/linklist.h \
				$(INCLUDE)/matrix.h $(INCLUDE)/BSQHydro.h $(SRC)/main.cpp
	$(CC) $(CFLAGS) $(SRC)/main.cpp

$(OBJDIR)/vector.o: $(INCLUDE)/vector.h $(INCLUDE)/mathdef.h $(SRC)/vector.cpp
	$(CC) $(CFLAGS) $(SRC)/vector.cpp

$(OBJDIR)/matrix.o: $(SRC)/matrix.cpp $(INCLUDE)/matrix.h $(INCLUDE)/vector.h $(INCLUDE)/mathdef.h
	$(CC) $(CFLAGS) $(SRC)/matrix.cpp

$(OBJDIR)/system_state.o: $(SRC)/system_state.cpp $(INCLUDE)/system_state.h
	$(CC) $(CFLAGS) $(SRC)/system_state.cpp

$(OBJDIR)/sph_workstation.o: $(SRC)/sph_workstation.cpp $(INCLUDE)/sph_workstation.h $(INCLUDE)/transport_coefficients.h
	$(CC) $(CFLAGS) $(SRC)/sph_workstation.cpp

$(OBJDIR)/linklist.o: $(SRC)/linklist.cpp $(INCLUDE)/linklist.h $(INCLUDE)/particle.h
	$(CC) $(CFLAGS) $(SRC)/linklist.cpp

$(OBJDIR)/input_output.o: $(SRC)/input_output.cpp $(INCLUDE)/input_output.h
	$(CC) $(CFLAGS) $(SRC)/input_output.cpp
	
$(OBJDIR)/kernel.o: $(SRC)/kernel.cpp $(INCLUDE)/kernel.h
	$(CC) $(CFLAGS) $(SRC)/kernel.cpp

$(OBJDIR)/BSQHydro.o: $(INCLUDE)/vector.h $(INCLUDE)/mathdef.h $(INCLUDE)/particle.h $(INCLUDE)/BSQHydro.h $(INCLUDE)/linklist.h $(INCLUDE)/eos.h $(INCLUDE)/matrix.h $(SRC)/BSQHydro.cpp
	$(CC) $(CFLAGS) $(SRC)/BSQHydro.cpp

$(OBJDIR)/particle.o: $(SRC)/particle.cpp $(INCLUDE)/particle.h $(INCLUDE)/eos.h
	$(CC) $(CFLAGS) $(SRC)/particle.cpp

$(OBJDIR)/eos.o: $(INCLUDE)/eos.h $(SRC)/eos.cpp $(INCLUDE)/eos_delaunay.h $(SRC)/eos_delaunay.cpp
	$(CC) $(CFLAGS) $(SRC)/eos.cpp

$(OBJDIR)/eos_table.o: $(INCLUDE)/eos_table.h $(SRC)/eos_table.cpp
	$(CC) $(CFLAGS) $(SRC)/eos_table.cpp

$(OBJDIR)/eos_extension.o: $(INCLUDE)/eos_extension.h $(SRC)/eos_extension.cpp
	$(CC) $(CFLAGS) $(SRC)/eos_extension.cpp

$(OBJDIR)/rootfinder.o: $(INCLUDE)/rootfinder.h $(SRC)/rootfinder.cpp
	$(CC) $(CFLAGS) $(SRC)/rootfinder.cpp

$(OBJDIR)/eos_delaunay.o: $(SRC)/eos_delaunay.cpp $(INCLUDE)/eos_delaunay.h
	$(CC) $(CFLAGS) $(INCFLAGS) -fpic -c $(SRC)/eos_delaunay.cpp 
	
$(OBJDIR)/transport_coefficients.o: $(SRC)/transport_coefficients.cpp $(INCLUDE)/transport_coefficients.h
	$(CC) $(CFLAGS) $(INCFLAGS) -fpic -c $(SRC)/transport_coefficients.cpp 
	
	


clean:
	-rm $(OBJDIR)/*.o persephone
