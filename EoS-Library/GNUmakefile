OBJS = test.o eos.o
CC = g++
DEBUG = -g
STANDARD= -std=c++11
WARNINGS= -w -Wall -Wextra -pedantic -Wno-conversion -Wcast-align -Wunused -Wshadow -Wpointer-arith
OPTIMIZATION= -Og
EIGENPATH=$(realpath ../eigen)
SPLINTERPATH=$(realpath ../splinter/include)
SPLINTERBUILDPATH=$(realpath ../splinter/build)
CFLAGS = $(WARNINGS) $(OPTIMIZATION) -c -I$(EIGENPATH) -I$(SPLINTERPATH) $(DEBUG) $(STANDARD)
LDFLAGS = -lm -lgsl -lblas -lm -lstdc++ -L$(SPLINTERBUILDPATH) -lsplinter-3-0 $(DEBUG) $(STANDARD)

test_EoS : $(OBJS)
	$(CC) -g $(OBJS) -o test_EoS $(LDFLAGS)

test.o: eos.h test.cpp eos.cpp
	$(CC) $(CFLAGS) test.cpp

eos.o: eos.h eos.cpp $(SPLINTERPATH)/datatable.h $(SPLINTERPATH)/bspline.h $(SPLINTERPATH)/bsplinebuilder.h
	$(CC) $(CFLAGS) eos.cpp

clean:
	-rm *.o test_EoS

