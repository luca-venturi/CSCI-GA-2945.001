/******************************************************************************
* FILE: omp_fixed4.c
* Luca Venturi comment: The omp_bug4 was not working due to a 'segmentation fault' (memory problem). I solved the problem by allocating the memory dinamically for each thread with 'malloc'.
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define N 1048

int main (int argc, char *argv[]) 
{
  
int nthreads, tid, i, j;
double *a;

/* Fork a team of threads with explicit variable scoping */
#pragma omp parallel shared(nthreads) private(i,j,tid,a)
  {
  /* Obtain/print thread info */
  tid = omp_get_thread_num();
  if (tid == 0) 
    {
    nthreads = omp_get_num_threads();
    printf("Number of threads = %d\n", nthreads);
    }
  printf("Thread %d starting...\n", tid);

  /* Each thread works on its own private copy of the array */
  a = (double *) malloc(N*N*sizeof(double));
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      a[i*N+j] = tid + i + j;

  /* For confirmation */
  printf("Thread %d done. Last element= %f\n",tid,a[N*N-1]);

  free(a);

  }  /* All threads join master thread and disband */

  return 0;
}

