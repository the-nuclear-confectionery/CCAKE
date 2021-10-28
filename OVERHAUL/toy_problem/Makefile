OBJS = main.o animal.o 
CC = g++
DEBUG = -g
STANDARD= -std=c++11
WARNINGS= -w -Wall -Wextra -pedantic -Wno-conversion -Wcast-align -Wunused -Wshadow -Wpointer-arith
OPTIMIZATION= -Og
CFLAGS = $(WARNINGS) $(OPTIMIZATION) -c - $(DEBUG) $(STANDARD)
LDFLAGS = -lstdc++ $(DEBUG) $(STANDARD)

animalCheck : $(OBJS)
	$(CC) -g $(OBJS) -o animalCheck $(LDFLAGS)

main.o: animal.h  main.cpp
	$(CC) $(CFLAGS) main.cpp

vector.o: animal.h animal.cpp
	$(CC) $(CFLAGS) animal.cpp
	
	

clean:
	-rm *.o animalCheck
