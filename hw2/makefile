CC = gcc
FLAGS = -O3 -fopenmp -Wall -g -lm
EXECS = omp_fixed2 omp_fixed3 omp_fixed4 omp_bug5 omp_bug6 jacobi2D

all: ${EXECS}

omp_fixed2: omp_fixed2.c
	${CC} ${FLAGS} $^ -o omp_fixed2

omp_fixed3: omp_fixed3.c
	${CC} ${FLAGS} $^ -o omp_fixed3

omp_fixed4: omp_fixed4.c
	${CC} ${FLAGS} $^ -o omp_fixed4

omp_bug5: omp_bug5.c
	${CC} ${FLAGS} $^ -o omp_bug5

omp_bug6: omp_bug6.c
	${CC} ${FLAGS} $^ -o omp_bug6

jacobi2D: jacobi2D.c
	${CC} ${FLAGS} $^ -o jacobi2D

clean:
	rm -f ${EXECS}
