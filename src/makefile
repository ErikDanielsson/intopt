CC=gcc
CFLAGS= -O3 -g -Werror -Wall
LDFLAGS= -lm -pthread
OBJ=bnb.o heap.o intopt.o main.o simplex.o worker.o

a.out: $(OBJ)
	$(CC) $(CFLAGS)  -o a.out $(OBJ) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o a.out
