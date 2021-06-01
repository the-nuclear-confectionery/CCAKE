# ==============================================================================
#  Purpose: Makefile for hydrodynamic interpolation and HDF converter
#  Author:  Christopher Plumberg
#  Date:    June 1, 2021
# ==============================================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			distclean	remove all binaries
##  

# HDF5 flags
export HDF5_CXX=$(which g++)
export HDF5_CXXLINKER=$(which g++)
export HDF5_CLINKER=$(which g++)

# Set compiler and flags
CC := h5c++
CFLAGS= -g -std=c++14 -lgsl -lgslcblas -lm
EIGENPATH=$(realpath ../eigen)
SPLINTERPATH=$(realpath ../splinter/include)
SPLINTERBUILDPATH=$(realpath ../splinter/build)
CFLAGS+= -c -I$(EIGENPATH) -I$(SPLINTERPATH)
LDFLAGS+= -lm -lgsl -lblas -lm -lstdc++ -L$(SPLINTERBUILDPATH) -lsplinter-3-0


# Various directories and definitions
RM          =   rm -f
O           =   .o
LDFLAGS     =   
INCFLAGS    =   
SYSTEMFILES =   $(SRCGNU)

# --------------- Files involved ------------------

MAIN		 =	test_EoS

MAINSRC      =  test.cpp

INC			 = 	eos.h
SRC			 =	eos.cpp

# -------------------------------------------------

TARGET		=	$(MAIN)

# --------------- Pattern rules -------------------

./%.o: ./%.cpp
	$(CC) $(CFLAGS) -c $< -o $@


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
eos.cpp:			  eos.h
test.cpp:             eos.h eos.cpp
