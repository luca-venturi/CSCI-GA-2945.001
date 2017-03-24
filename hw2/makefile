CC = gcc
FLAGS = -O3 -fopenmp -Wall -g
EXECS = omp_solved2 omp_solved3 omp_solved4 omp_solved5 omp_solved6 jacobi2D gs2D jacobi2D-omp gs2D-omp

all: ${EXECS}

omp_solved2: omp_solved2.c
	${CC} ${FLAGS} $^ -o omp_solved2

omp_solved3: omp_solved3.c
	${CC} ${FLAGS} $^ -o omp_solved3

omp_solved4: omp_solved4.c
	${CC} ${FLAGS} $^ -o omp_solved4

omp_solved5: omp_solved5.c
	${CC} ${FLAGS} $^ -o omp_solved5

omp_solved6: omp_solved6.c
	${CC} ${FLAGS} $^ -o omp_solved6

jacobi2D: jacobi2D.c
	${CC} ${FLAGS} $^ -o jacobi2D -lm

jacobi2D-omp: jacobi2D-omp.c
	${CC} ${FLAGS} $^ -o jacobi2D-omp -lm

gs2D: gs2D.c
	${CC} ${FLAGS} $^ -o gs2D -lm

gs2D-omp: gs2D-omp.c
	${CC} ${FLAGS} $^ -o gs2D-omp -lm

clean:
	rm -f ${EXECS}
