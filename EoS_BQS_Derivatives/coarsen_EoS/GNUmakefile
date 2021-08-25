# Set compiler and flags
CC := g++
CFLAGS= -O3 -std=c++11 -lgsl -lgslcblas -lm -g

# Various directories and definitions
RM          =   rm -f
O           =   .o
LDFLAGS     =   -L/usr/local/src/gsl/2.5/lib -L/projects/jnorhos/plumberg/qhull/lib
INCFLAGS    =   -I/usr/local/src/gsl/2.5/include -I/projects/jnorhos/BSQ/kdtree_test \
				-I/projects/jnorhos/plumberg/qhull/src
SYSTEMFILES =   $(SRCGNU)


# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	coarsen_EoS
endif
ifeq "$(MAIN2)" ""
MAIN2		=	coarsen_EoS_faster
endif

MAINSRC     =   coarsen_EoS.cpp
INC		= 	
MAIN2SRC     =   coarsen_EoS_faster.cpp
INC2		= 	Stopwatch.h

# -------------------------------------------------

TARGET		=	$(MAIN)
TARGET2		=	$(MAIN2)
INSTPATH	=	..

# --------------- Pattern rules -------------------

$(TARGET):
	$(CC) $(MAINSRC) -o $(TARGET) $(CFLAGS) $(INCFLAGS) $(LDFLAGS)

$(TARGET2):
	$(CC) $(MAIN2SRC) -o $(TARGET2) $(CFLAGS) $(INCFLAGS) $(LDFLAGS)

interpolate_ebsq_v2: interpolate_ebsq_v2.cpp
	$(CC) interpolate_ebsq_v2.cpp -o interpolate_ebsq_v2 $(CFLAGS) $(INCFLAGS) -fpic \
			$(LDFLAGS) -lqhullcpp -lqhullstatic_r

eos_delaunay.o: eos_delaunay.cpp eos_delaunay.h
	$(CC) $(CFLAGS) -c eos_delaunay.cpp $(INCFLAGS) -fpic \
			$(LDFLAGS) -lqhullcpp -lqhullstatic_r

interpolate_ebsq: eos_delaunay.o
	$(CC) -g eos_delaunay.o interpolate_ebsq.cpp -o interpolate_ebsq $(CFLAGS) \
			$(INCFLAGS) -fpic $(LDFLAGS) -lqhullcpp -lqhullstatic_r

# -------------------------------------------------

.PHONY:		all help distclean install

all:		$(TARGET) $(TARGET2) interpolate_ebsq_v2 interpolate_ebsq

help:
		@grep '^##' GNUmakefile

distclean:	
		-rm $(TARGET) $(TARGET2)

install:	$(TARGET) $(TARGET2)
		cp $(TARGET) $(INSTPATH)
		cp $(TARGET2) $(INSTPATH)

# --------------- Dependencies -------------------
coarsen_EoS.cpp:		
coarsen_EoS_faster.cpp:		Stopwatch.h
interpolate_ebsq_v2.cpp:	Stopwatch.h delaunay.h point_in_simplex.h
eos_delaunay.cpp:			eos_delaunay.h
interpolate_ebsq.cpp:		eos_delaunay.cpp eos_delaunay.h
