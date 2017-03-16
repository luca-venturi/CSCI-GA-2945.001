/******************************************************************************
* FILE: omp_fixed6.c
* Luca Venturi comment: The bug was caused by a wrong declaration of the variable 'sum'. When the function was called by every thread a private variable 'sum' was created in every thread: being private, this variable could not be declared as reduction. I adjusted the program by changing the code so that 'sum' is actually a reduction sum.
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define VECLEN 100

float a[VECLEN], b[VECLEN];

int main (int argc, char *argv[]) {
int i, tid;
float sum=0;

for (i=0; i < VECLEN; i++)
  a[i] = b[i] = 1.0 * i;

#pragma omp parallel private(tid) reduction(+:sum)
  {
  tid = omp_get_thread_num();

  #pragma omp for
    for (i=0; i < VECLEN; i++) {
      sum = sum + (a[i]*b[i]);
      printf("  tid=%d i=%d\n",tid,i);
    }
  }

printf("Sum = %f\n",sum);

return 0;
}
