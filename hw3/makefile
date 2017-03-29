CC = mpicc
FLAGS = -O3 -Wall -g
EXECS = int_ring int_ring_largearray

all: ${EXECS}

int_ring: int_ring.c
	${CC} ${FLAGS} $^ -o int_ring

int_ring_largearray: int_ring_largearray.c
	${CC} ${FLAGS} $^ -o int_ring_largearray

clean:
	rm -f ${EXECS}
