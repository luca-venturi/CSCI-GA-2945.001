EXECS = multigrid
CC = gcc
FLAGS = -O3 -Wall -g -lm

all: ${EXECS}

multigrid: multigrid.c
	${CC} ${FLAGS} $^ -o multigrid

multigrid-mpi: multigrid-mpi.c
	mpicc ${FLAGS} $^ -o multigrid-mpi

clean:
	rm -f ${EXECS}
