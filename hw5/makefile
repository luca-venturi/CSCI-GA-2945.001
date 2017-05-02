EXECS = multigrid
CC = gcc
FLAGS = -O3 -Wall -g -lm

all: ${EXECS}

multigrid: multigrid.c
	${CC} ${FLAGS} $^ -o multigrid

clean:
	rm -f ${EXECS}
