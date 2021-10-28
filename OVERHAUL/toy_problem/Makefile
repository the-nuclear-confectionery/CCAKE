OBJS = main.o animal.o 
#CC = g++
CC = h5c++
DEBUG = -g
STANDARD= -std=c++11
WARNINGS= -w -Wall -Wextra -pedantic -Wno-conversion -Wcast-align -Wunused -Wshadow -Wpointer-arith
OPTIMIZATION= -Og
#OPTIMIZATION= -Ofast -funroll-loops -fvectorize -finline-functions -s -DNDEBUG -march=native -mtune=native
# change to -O3 or -Ofast for production
CFLAGS = $(WARNINGS) $(OPTIMIZATION) -c -I$(EOSLIBPATH) -I$(EIGENPATH) \
        -I$(GSLPATH) $(DEBUG) $(STANDARD)
LDFLAGS = -L$(EOSLIBPATH) -lbsqeos -lm -lgsl -lgslcblas -lblas -lstdc++ \
          -L$(GSLLIBPATH) $(DEBUG) $(STANDARD)

export HDF5_LDFLAGS+=$(LDFLAGS)

animalCheck : $(OBJS)
	$(CC) -g $(OBJS) -o animalChecl $(LDFLAGS)

main.o: animal.h  main.cpp
	$(CC) $(CFLAGS) main.cpp

vector.o: animal.h animal.cpp
	$(CC) $(CFLAGS) animal.cpp
	
	

clean:
	-rm *.o animalCheck
