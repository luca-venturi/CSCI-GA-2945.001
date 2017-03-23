/******************************************************************************
* FILE: omp_fixed6.c
* Luca Venturi comment: The bug was caused by a wrong declaration of the variable 'sum'. The problem coulb fixed in many different ways. In my solved version, all the parallel pragma is in the function 'dotprod()', and sum is a reduction variable there, which is returned at the end by the function to the 'main()' program.
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define VECLEN 100

float a[VECLEN], b[VECLEN];

float dotprod ()
{
int i,tid;
float sum=0.0;

#pragma omp parallel default(none) private(i,tid) shared(a,b) reduction(+:sum)
  {
tid = omp_get_thread_num();
#pragma omp for 
  for (i=0; i < VECLEN; i++)
    {
    sum = sum + (a[i]*b[i]);
    printf("  tid= %d i=%d\n",tid,i);
    }
  }

  return sum;
}


int main (int argc, char *argv[]) {
int i;
float sum;

for (i=0; i < VECLEN; i++)
  a[i] = b[i] = 1.0 * i;

sum = dotprod();

printf("Sum = %f\n",sum);

return 0;
}

