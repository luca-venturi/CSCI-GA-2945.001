CC = gcc
FLAGS = -O3 -fopenmp -Wall -g
EXECS = omp_fixed2 omp_fixed3 omp_fixed4 omp_fixed5 omp_fixed6 jacobi2D gs2D jacobi2D-omp gs2D-omp

all: ${EXECS}

omp_fixed2: omp_fixed2.c
	${CC} ${FLAGS} $^ -o omp_fixed2

omp_fixed3: omp_fixed3.c
	${CC} ${FLAGS} $^ -o omp_fixed3

omp_fixed4: omp_fixed4.c
	${CC} ${FLAGS} $^ -o omp_fixed4

omp_fixed5: omp_fixed5.c
	${CC} ${FLAGS} $^ -o omp_fixed5

omp_fixed6: omp_fixed6.c
	${CC} ${FLAGS} $^ -o omp_fixed6

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
