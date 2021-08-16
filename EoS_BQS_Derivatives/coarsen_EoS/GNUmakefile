# Set compiler and flags
CC := g++
CFLAGS= -std=c++11 -lgsl -lgslcblas -lm -g

# Various directories and definitions
RM          =   rm -f
O           =   .o
LDFLAGS     =   -L/usr/local/src/gsl/2.5/lib
INCFLAGS    =   -I/usr/local/src/gsl/2.5/include
SYSTEMFILES =   $(SRCGNU)


# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	coarsen_EoS
endif

MAINSRC     =   coarsen_EoS.cpp

INC		= 	

# -------------------------------------------------

TARGET		=	$(MAIN)
INSTPATH	=	..

# --------------- Pattern rules -------------------

$(TARGET):
	$(CC) $(MAINSRC) -o $(TARGET) $(CFLAGS) $(INCFLAGS)  $(LDFLAGS)

# -------------------------------------------------

.PHONY:		all help distclean install

all:		$(TARGET)

help:
		@grep '^##' GNUmakefile

distclean:	
		-rm $(TARGET)

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
coarsen_EoS.cpp:		

