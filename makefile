CC = gcc
FLAGS = -O3 -fopenmp -Wall -g
EXECS = omp_bug2 omp_bug3 omp_bug4 omp_bug5 omp_bug6

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
