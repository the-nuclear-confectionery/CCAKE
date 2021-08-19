# Set compiler and flags
CC := g++
CFLAGS= -O3 -std=c++11 -lgsl -lgslcblas -lm -g

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
	$(CC) $(MAINSRC) -o $(TARGET) $(CFLAGS) $(INCFLAGS)  $(LDFLAGS)

$(TARGET2):
	$(CC) $(MAIN2SRC) -o $(TARGET2) $(CFLAGS) $(INCFLAGS)  $(LDFLAGS)

interpolate_ebsq: interpolate_ebsq.cpp Stopwatch.h
	$(CC) interpolate_ebsq.cpp -o $(TARGET2) $(CFLAGS) $(INCFLAGS)  $(LDFLAGS)

# -------------------------------------------------

.PHONY:		all help distclean install

all:		$(TARGET) $(TARGET2) interpolate_ebsq

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
interpolate_ebsq.cpp:		Stopwatch.h

