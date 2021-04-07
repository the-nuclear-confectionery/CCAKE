OBJS = main.o vector.o hydrosim.o newspline.o random.o eos.o matrix.o enteric.o
CC = g++
DEBUG = -g -std=c++11 
EIGENPATH=$(realpath ../eigen)
SPLINTERPATH=$(realpath ../splinter/include)
SPLINTERBUILDPATH=$(realpath ../splinter/build)
CFLAGS = -Wall -c -I$(EIGENPATH) -I$(SPLINTERPATH) $(DEBUG)
LDFLAGS = -lm -O3 -lgsl -lblas -lm -lstdc++ -L$(SPLINTERBUILDPATH) -lsplinter-3-0 $(DEBUG)

vusphydro : $(OBJS)
	$(CC) -g $(OBJS) -o vusphydro $(LDFLAGS)

main.o: vector.h mathdef.h particle.h eos.h enteric.h LinkList.h random.h matrix.h hydrosim.h eostables.h main.cpp
	$(CC) $(CFLAGS) main.cpp

vector.o: vector.h mathdef.h vector.cpp
	$(CC) $(CFLAGS) vector.cpp

matrix.o: matrix.cpp matrix.h vector.h mathdef.h
	$(CC) $(CFLAGS) matrix.cpp

hydrosim.o: vector.h mathdef.h particle.h hydrosim.h LinkList.h rungekutta4.h eos.h output.h bbmg.h matrix.h hydrosim.cpp
	$(CC) $(CFLAGS) hydrosim.cpp

newspline.o: newspline.cpp newspline.hpp
	$(CC) $(CFLAGS) newspline.cpp

eos.o: tables.h eos.h eostables.h eos.cpp $(SPLINTERPATH)/datatable.h $(SPLINTERPATH)/bspline.h $(SPLINTERPATH)/bsplinebuilder.h
	$(CC) $(CFLAGS) eos.cpp

random.o: random.cpp random.h
	$(CC) $(CFLAGS) random.cpp

enteric.o:  vector.h particle.h eos.h tables.h enteric.h LinkList.h matrix.h eostables.h enteric.cpp
	$(CC) $(CFLAGS) enteric.cpp



# deletes all .o files in case of error in makefile.  This should not be something needed regularly, only if you are updating the code and run into an error.
#clean:
#	\rm *.o *~ vusphydro

# puts all the files in a vusphydro.tar file
#tar:
#    tar cfv vusphydro.tar vector.h main.cpp vector.cpp
