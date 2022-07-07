# Set compiler and flags
CC := h5c++
CFLAGS= -O3 -std=c++11 -g

# Various directories and definitions
RM          =   rm -f
O           =   .o
LDFLAGS     =   -L/projects/jnorhos/plumberg/qhull/lib -lqhullcpp -lqhullstatic_r \
		-L/projects/jnorhos/local/lib -lgsl -lgslcblas -lm \
		-L`realpath ../../EoS_BQS_Derivatives/lib` -lbsqeos
INCFLAGS    =   -I/projects/jnorhos/plumberg/qhull/src \
		-I/projects/jnorhos/local/include \
		-I`realpath ../../EoS_BQS_Derivatives/lib`
SYSTEMFILES =   $(SRCGNU)
OBJS	    =	eos_delaunay.o interpolate_ebsq.o delaunay.o point_in_simplex.o

# -------------------------------------------------

TARGET		=	interpolate_ebsq
INSTPATH	=	..

# --------------- Pattern rules -------------------

$(TARGET): $(OBJS)
	$(CC) -g $(OBJS) -o interpolate_ebsq $(LDFLAGS)

eos_delaunay.o: eos_delaunay.cpp eos_delaunay.h
	$(CC) $(CFLAGS) $(INCFLAGS) -fpic -c eos_delaunay.cpp 

delaunay.o: delaunay.cpp delaunay.h
	$(CC) $(CFLAGS) -c delaunay.cpp $(INCFLAGS) -fpic

point_in_simplex.o: point_in_simplex.cpp point_in_simplex.h
	$(CC) $(CFLAGS) -c point_in_simplex.cpp $(INCFLAGS) -fpic

interpolate_ebsq.o: interpolate_ebsq.cpp eos_delaunay.h
	$(CC) $(CFLAGS) -c interpolate_ebsq.cpp $(INCFLAGS) -fpic


# -------------------------------------------------

.PHONY:		all help clean install

all:		$(TARGET)

help:
		@grep '^##' GNUmakefile

clean:	
		-rm $(TARGET) $(OBJS)

install: $(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
eos_delaunay.cpp:			eos_delaunay.h
interpolate_ebsq.cpp:		eos_delaunay.cpp eos_delaunay.h Stopwatch.h delaunay.h point_in_simplex.h
