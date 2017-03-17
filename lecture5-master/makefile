CC = gcc
FLAGS = -O3 -fopenmp -Wall -g
EXECS = jacobi jacobi-omp

all: ${EXECS}

jacobi: jacobi.c
	${CC} ${FLAGS} $^ -o jacobi -lm

jacobi-omp: jacobi-omp.c
	${CC} ${FLAGS} $^ -o jacobi-omp -lm

gs2D: gs2D.c
	${CC} ${FLAGS} $^ -o gs2D -lm

clean:
	rm -f ${EXECS}
