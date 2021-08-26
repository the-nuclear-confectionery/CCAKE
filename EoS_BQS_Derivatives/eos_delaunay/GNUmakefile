# Set compiler and flags
CC := g++
CFLAGS= -O3 -std=c++11

# Various directories and definitions
RM          =   rm -f
O           =   .o
LDFLAGS     =   -L/usr/local/src/gsl/2.5/lib -L/projects/jnorhos/plumberg/qhull/lib \
				-lgsl -lgslcblas -lm
INCFLAGS    =   -I/usr/local/src/gsl/2.5/include \
				-I/projects/jnorhos/plumberg/qhull/src
SYSTEMFILES =   $(SRCGNU)
OBJS		=	eos_delaunay.o interpolate_ebsq.o

# -------------------------------------------------

TARGET		=	interpolate_ebsq
INSTPATH	=	..

# --------------- Pattern rules -------------------

$(TARGET): $(OBJS)
	$(CC) -g $(OBJS) -o interpolate_ebsq -fpic $(LDFLAGS) -lqhullcpp -lqhullstatic_r

eos_delaunay.o: eos_delaunay.cpp eos_delaunay.h
	$(CC) $(CFLAGS) -c eos_delaunay.cpp $(INCFLAGS) -fpic

interpolate_ebsq.o: interpolate_ebsq.cpp eos_delaunay.h
	$(CC) $(CFLAGS) -c interpolate_ebsq.cpp $(INCFLAGS) -fpic


# -------------------------------------------------

.PHONY:		all help distclean install

all:		$(TARGET)

help:
		@grep '^##' GNUmakefile

distclean:	
		-rm $(TARGET)

install: $(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
eos_delaunay.cpp:			eos_delaunay.h
interpolate_ebsq.cpp:		eos_delaunay.cpp eos_delaunay.h Stopwatch.h delaunay.h point_in_simplex.h
