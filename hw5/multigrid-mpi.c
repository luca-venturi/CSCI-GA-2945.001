/* */
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "util.h"
#include <string.h>

/* compute norm of global vector */
double compute_norm(double *lu, int lN)
{
	int i;
	double gres = 0.0, lres = 0.0;

	for (i = 1; i <= lN + 1; i++)
		lres += lu[i] * lu[i];
  
	/* use allreduce for convenience; a reduce would also be sufficient */
	MPI_Allreduce(&lres, &gres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return sqrt(gres);
}

/* set vector to zero */
void set_zero (double *lu, int lN) {
	int i;
	for (i = 0; i <= lN + 1; i++)
		lu[i] = 0.0;
}

/* debug function */
void output_to_screen (double *lu, int lN) {
	int i, j;
	if (0 == rank)
		printf("%f ", lu[0]);
	for (j = 0; j < size; j++) {
		if (j == rank) {
			for (i = 1; i <= lN; i++)
				printf("%f ", lu[i]);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	if (size - 1 == rank)
		printf("%f ", lu[lN+1]);
	printf("\n");
}

/* coarsen uf from length lN+2 to lenght lN/2+2
   assuming N = 2^l
*/
void coarsen(double *luf, double *luc, int lN) {
  int ic;
  for (ic = 1; ic < N/2; ++ic)
    uc[ic] = 0.5 * uf[2*ic] + 0.25 * (uf[2*ic-1]+uf[2*ic+1]);
}

