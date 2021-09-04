TARGET = libbsqeos.a

CC=h5cc
FLAGS_NEEDED_FOR_GCC_10_3_0 = -fcommon
INCFLAGS=-I/projects/jnorhos/local/include
LDFLAGS=$(FLAGS_NEEDED_FOR_GCC_10_3_0) -L/projects/jnorhos/local/lib -lm -lgsl -lgslcblas -lblas
#SRC=../src/lnsrchD.c ../src/fminD.c ../src/fdjacD.c ../src/newtD.c ../src/Functions.c ../src/lubksbD.c ../src/ludcmpD.c ../src/nrutilD.c ../src/nm_solver.c

OBJS=driver.o ../src/lnsrchD.o ../src/fminD.o ../src/fdjacD.o ../src/newtD.o ../src/Functions.o ../src/lubksbD.o ../src/ludcmpD.o ../src/nrutilD.o ../src/nm_solver.o

$(TARGET): $(OBJS)
	ar rcs $@ $^

main: $(TARGET) main.o
	$(CC) $^ -o $@ $(LDFLAGS)

$(OBJS): %.o: %.c
	$(CC) -c -o $@ $< $(INCFLAGS) $(LDFLAGS)

clean:
	rm -f *.o $(TARGET)
