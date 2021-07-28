# Set compiler and flags
CC := g++
GSLPATH=/projects/jnorhos/local/gsl/include
GSLLIBPATH=/projects/jnorhos/local/gsl/lib
CFLAGS= -g -std=c++14 -I$(GSLPATH) -L$(GSLLIBPATH) -lm -lgsl -lgslcblas -lblas

# Various directories and definitions
RM          =   rm -f
O           =   .o
LDFLAGS     =   
INCFLAGS    =   
SYSTEMFILES =   $(SRCGNU)

# --------------- Files involved ------------------

MAIN		 =	interp_thermo

MAINSRC      =  interp_thermo.cpp

INC			 = 	interp_thermo.h Stopwatch.h

# -------------------------------------------------

TARGET		=	$(MAIN)

# --------------- Pattern rules -------------------

$(TARGET):
	$(CC) $(MAINSRC) -o $(TARGET) $(CFLAGS) $(INCFLAGS) $(LDFLAGS)

# -------------------------------------------------

.PHONY:		all help distclean

all:		$(TARGET)

help:
		@grep '^##' GNUmakefile

distclean:	
		-rm $(TARGET)

# --------------- Dependencies -------------------
interp_thermo.cpp:	interp_thermo.h Stopwatch.h
