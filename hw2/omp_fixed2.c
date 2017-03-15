/******************************************************************************
* FILE: omp_fixed2.c
* DESCRIPTION: The bug was caused because 'total' and 'tid' variables were considered as shared variables by default (since nothing was specified). Clearly 'tid' has to be set as private. I also set 'total' as private: in this way every thread computes a part of the total sum N*(N-1)/2. Another option would be to set . Also I changed the type of 'total' from 'float' to 'long double' (in this way the result we get is correct). 
* AUTHOR: Luca Venturi
******************************************************************************/ 
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) 
{
int nthreads, i, tid;
long double total;

/*** Spawn parallel region ***/
#pragma omp parallel private(tid, total)
  {
  /* Obtain thread number */
  tid = omp_get_thread_num();
  /* Only master thread does this */
  if (tid == 0) {
    nthreads = omp_get_num_threads();
    printf("Number of threads = %d\n", nthreads);
    }
  printf("Thread %d is starting...\n",tid);

  #pragma omp barrier

  /* do some work */
  total=0.0;
#pragma omp for schedule(dynamic,10)
  for (i=0; i<1000000; i++) 
     total +=  i*1.0;

  printf("Thread %d is done! Total= %Le\n",tid,total);

  } /*** End of parallel region ***/

  return 0;
}
