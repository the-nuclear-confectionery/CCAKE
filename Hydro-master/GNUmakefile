# ==============================================================================
#  Purpose: Makefile for hydrodynamic interpolation and HDF converter
#  Author:  Christopher Plumberg
#  Date:    February 2, 2021
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

# Various directories and definitions
RM          =   rm -f
O           =   .o
LDFLAGS     =   
INCFLAGS    =   
SYSTEMFILES =   $(SRCGNU)

# --------------- Files involved ------------------

MAIN		 =	convert_to_hdf

MAINSRC      =  convert_to_hdf.cpp

INC			 = 	convert_to_hdf.h

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
convert_to_hdf.cpp:             convert_to_hdf.h
