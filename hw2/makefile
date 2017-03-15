CC = gcc
FLAGS = -O3 -fopenmp -Wall -g -lm
EXECS = omp_bug2 omp_bug3 omp_bug4 omp_bug5 omp_bug6 jacobi2D

all: ${EXECS}

omp_bug2: omp_bug2.c
	${CC} ${FLAGS} $^ -o omp_bug2

omp_bug3: omp_bug3.c
	${CC} ${FLAGS} $^ -o omp_bug3

omp_bug4: omp_bug4.c
	${CC} ${FLAGS} $^ -o omp_bug4

omp_bug5: omp_bug5.c
	${CC} ${FLAGS} $^ -o omp_bug5

omp_bug6: omp_bug6.c
	${CC} ${FLAGS} $^ -o omp_bug6

jacobi2D: jacobi2D.c
	${CC} ${FLAGS} $^ -o jacobi2D
