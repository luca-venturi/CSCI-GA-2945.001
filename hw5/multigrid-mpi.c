/* */
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "util.h"
#include <string.h>

/* returns 1 if N = 2^k for some k=0,1,..., 0 otherwise */
int pow_of_two(int N)
{
	if ((N & (N-1)) == 0)
		return 1;
	else
		return 0;
}

/* compute norm of global vector */
double compute_norm(double *lu, int lN)
{
	int i;
	double gres = 0.0, lres = 0.0;

	for (i = 1; i <= lN; i++)
		lres += lu[i] * lu[i];
  
	/* use allreduce for convenience; a reduce would also be sufficient */
	MPI_Allreduce(&lres, &gres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return sqrt(gres);
}

/* set vector to zero */
void set_zero (double *lu, int lN) {
	int i;
	for (i = 0; i <= lN + 3; i++)
		lu[i] = 0.0;
}

/* debug function */
void output_to_screen (double *lu, int lN, int rank, int size) {
	int i, j;
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

/* coarsen uf from length lN+3 to lenght lN/2+3 assuming N = 2^l */
void coarsen(double *luf, double *luc, int lN, int rank, int size) {
	int ic, lNhalf = lN/2;
	MPI_Request request[4];
	for (ic = 1; ic <= lNhalf+1; ++ic)
		luc[ic] = 0.5 * luf[2*ic-1] + 0.25 * (luf[2*ic-2]+luf[2*ic]);
	/* communicate ghost values */
	if (rank < size - 1) {
      MPI_Isend(&luc[lNhalf], 1, MPI_DOUBLE, rank+1, 124, MPI_COMM_WORLD, request);
      MPI_Irecv(&luc[lNhalf+2], 1, MPI_DOUBLE, rank+1, 123, MPI_COMM_WORLD, request+1);
  }
	if (rank > 0) {
      MPI_Isend(&luc[2], 1, MPI_DOUBLE, rank-1, 123, MPI_COMM_WORLD, request+2);
      MPI_Irecv(&luc[0], 1, MPI_DOUBLE, rank-1, 124, MPI_COMM_WORLD, request+3);
  }
	MPI_Barrier(MPI_COMM_WORLD);
	/* set boundary values to 0 */
	if (rank == 0)		
		luc[0] = luc[1] = 0.0;
	else if (rank == size-1)		
		luc[lNhalf+1] = luc[lNhalf+2] = 0.0;
}

/* refine lu from length lN+3 to lenght 2*lN+3 assuming N = 2^l, and add to existing luf; ghost values do not change */
void refine_and_add(double *lu, double *luf, int lN, int rank, int size)
{
  int i;
	MPI_Request request[4];
	/* update values */
  for (i = 1; i <= lN; ++i) {
    luf[2*i-1] += lu[i];
		luf[2*i] += 0.5 * (lu[i] + lu[i+1]);
  }
	luf[2*lN+1] += lu[lN+1];
	/* update ghost values */
	if (rank < size - 1) {
      MPI_Isend(&luf[2*lN], 1, MPI_DOUBLE, rank+1, 124, MPI_COMM_WORLD, request);
      MPI_Irecv(&luf[2*lN+2], 1, MPI_DOUBLE, rank+1, 123, MPI_COMM_WORLD, request+1);
  }
	if (rank > 0) {
      MPI_Isend(&luf[2], 1, MPI_DOUBLE, rank-1, 123, MPI_COMM_WORLD, request+2);
      MPI_Irecv(&luf[0], 1, MPI_DOUBLE, rank-1, 124, MPI_COMM_WORLD, request+3);
  }
	MPI_Barrier(MPI_COMM_WORLD);
	/* set boundary values to 0 */
	if (rank == 0)		
		luf[0] = luf[1] = 0.0;
	else if (rank == size-1)		
		luf[lN+1] = luf[lN+2] = 0.0;
}

/* compute residual vector */
void compute_residual(double *lu, double *lrhs, double *lres, int lN, double invhsq, int rank, int size)
{
	int i;
	MPI_Request request[4];
	/* update values */
	for (i = 1; i <= lN+1; i++)
		lres[i] = (lrhs[i] - (2.*lu[i] - lu[i-1] - lu[i+1]) * invhsq);
	/* communicate ghost values */
	if (rank < size - 1) {
      MPI_Isend(&lres[lN], 1, MPI_DOUBLE, rank+1, 124, MPI_COMM_WORLD, request);
      MPI_Irecv(&lres[lN+2], 1, MPI_DOUBLE, rank+1, 123, MPI_COMM_WORLD, request+1);
    }
	if (rank > 0) {
      MPI_Isend(&lres[2], 1, MPI_DOUBLE, rank-1, 123, MPI_COMM_WORLD, request+2);
      MPI_Irecv(&lres[0], 1, MPI_DOUBLE, rank-1, 124, MPI_COMM_WORLD, request+3);
    }
	MPI_Barrier(MPI_COMM_WORLD);
}

/* compute residual and coarsen */
void compute_and_coarsen_residual(double *lu, double *lrhs, double *lresc, int lN, double invhsq, int rank, int size)
{
	double *lresf = calloc(sizeof(double), lN+3);
	compute_residual(lu, lrhs, lresf, lN, invhsq, rank, size);
	coarsen(lresf, lresc, lN, rank, size);
	free(lresf);
}

/* Perform Jacobi iterations on u */
void jacobi(double *lu, double *lrhs, int lN, double hsq, int ssteps, int rank, int size)
{
	int i, j;
	MPI_Request request[4];
	/* Jacobi damping parameter -- plays an important role in MG */
	double omega = 2./3.;
	double *lunew = calloc(sizeof(double), lN+3);
	for (j = 0; j < ssteps; ++j) {
		/* update values */
		for (i = 1; i <= lN+1; i++)
			lunew[i]  = lu[i] +  omega * 0.5 * (hsq*lrhs[i] + lu[i-1] + lu[i+1] - 2*lu[i]);
		/* communicate ghost values */
		if (rank < size - 1) {
			MPI_Isend(&lunew[lN], 1, MPI_DOUBLE, rank+1, 124, MPI_COMM_WORLD, request);
			MPI_Irecv(&lunew[lN+2], 1, MPI_DOUBLE, rank+1, 123, MPI_COMM_WORLD, request+1);
		}
		if (rank > 0) {
			MPI_Isend(&lunew[2], 1, MPI_DOUBLE, rank-1, 123, MPI_COMM_WORLD, request+2);
			MPI_Irecv(&lunew[0], 1, MPI_DOUBLE, rank-1, 124, MPI_COMM_WORLD, request+3);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		/* set boundary values to 0 */
		if (rank == 0)		
			lunew[0] = lunew[1] = 0.0;
		else if (rank == size-1)		
			lunew[lN+1] = lunew[lN+2] = 0.0;
		/* memcpy */
		memcpy(lu, lunew, (lN+3)*sizeof(double));
	}
	free(lunew);
}

int main(int argc, char * argv[])
{
	int i, rank, size, lNfine, l, iter, max_iters, levels, ssteps = 3;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int size_err = pow_of_two(size);
	if (rank == 0 && !(size_err)) {
		fprintf(stderr, "size: # of processors must be power of two number\n");
		abort();
	}	

	/* get name of host running MPI process */
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(processor_name, &name_len);
	printf("Rank %d/%d running on %s.\n", rank, size, processor_name);

	if (argc < 3 || argc > 4) {
		fprintf(stderr, "Usage: ./multigrid_1d Nfine maxiter [s-steps]\n");
		fprintf(stderr, "Nfine: # of intervals, must be power of two number\n");
		fprintf(stderr, "s-steps: # jacobi smoothing steps (optional, default is 3)\n");
		abort();
	}
	sscanf(argv[1], "%d", &lNfine);
	sscanf(argv[2], "%d", &max_iters);
	if (argc > 3)
		sscanf(argv[3], "%d", &ssteps);

	/* compute number of multigrid levels */
	levels = floor(log2(lNfine));
	printf("Multigrid Solve using V-cycles for -u'' = f on (0,1)\n");
	printf("Number of intervals = %d, max_iters = %d\n", lNfine, max_iters);
	printf("Number of MG levels: %d \n", levels);

	/* timing */
	MPI_Barrier(MPI_COMM_WORLD);
	timestamp_type time1, time2;
	get_timestamp(&time1);

	/* Allocation of vectors, including left and right bdry points */
	double *lu[levels], *lrhs[levels];
	/* N, h*h and 1/(h*h) on each level */
	int *lN = (int*) calloc(sizeof(int), levels);
	double *invhsq = (double* ) calloc(sizeof(double), levels);
	double *hsq = (double* ) calloc(sizeof(double), levels);
	double *lres = (double *) calloc(sizeof(double), lNfine+3);
	for (l = 0; l < levels; ++l) {
		lN[l] = lNfine / (int) pow(2,l);
		double h = 1.0 / (size * lN[l]);
		hsq[l] = h * h;
		printf("MG level %2d, N = %8d\n", l, lN[l]);
		invhsq[l] = 1.0 / hsq[l];
		lu[l] = (double *) calloc(sizeof(double), lN[l]+3);
		lrhs[l] = (double *) calloc(sizeof(double), lN[l]+3);
	}
	/* rhs on finest mesh */
	for (i = 0; i <= lN[0]; ++i) {
		lrhs[0][i] = 1.0;
	}
	/* set boundary values (unnecessary if calloc is used) */
	if (rank == 0)
		lu[0][0] = lu[0][1] = 0.0;
	else if (rank == size-1)
		lu[0][lN[0]+1] = lu[0][lN[0]+2] = 0.0;
	double res_norm, res0_norm, tol = 1e-6;

	/* initial residual norm */
	compute_residual(lu[0], lrhs[0], lres, lN[0], invhsq[0], rank, size);
	res_norm = res0_norm = compute_norm(lres, lN[0]);
	printf("Initial Residual: %f\n", res0_norm); 

	for (iter = 0; iter < max_iters && res_norm/res0_norm > tol; iter++) {
		/* V-cycle: Coarsening */
		for (l = 0; l < levels-1; ++l) {
			/* pre-smoothing and coarsen */
			jacobi(lu[l], lrhs[l], lN[l], hsq[l], ssteps, rank, size);
			compute_and_coarsen_residual(lu[l], lrhs[l], lrhs[l+1], lN[l], invhsq[l], rank, size);
			/* initialize correction for solution with zero */
			set_zero(lu[l+1],lN[l+1]);
		}
		/* V-cycle: Solve on coarsest grid using many smoothing steps */
		jacobi(lu[levels-1], lrhs[levels-1], lN[levels-1], hsq[levels-1], 50, rank, size);

		/* V-cycle: Refine and correct */
		for (l = levels-1; l > 0; --l) {
			/* refine and add to u */
			refine_and_add(lu[l], lu[l-1], lN[l], rank, size);
			/* post-smoothing steps */
			jacobi(lu[l-1], lrhs[l-1], lN[l-1], hsq[l-1], ssteps, rank, size);
		}

		if (0 == (iter % 1)) {
			compute_residual(lu[0], lrhs[0], lres, lN[0], invhsq[0], rank, size);
			res_norm = compute_norm(lres, lN[0]);
			if (0 == rank)
				printf("[Iter %d] Residual norm: %2.8f\n", iter, res_norm);
		}
	}

	/* Clean up */
	free(hsq);
	free(invhsq);
	free(lN);
	free(lres);
	for (l = levels-1; l >= 0; --l) {
		free(lu[l]);
		free(lrhs[l]);
	}

	/* timing */
	MPI_Barrier(MPI_COMM_WORLD);
	get_timestamp(&time2);
	double elapsed = timestamp_diff_in_seconds(time1,time2);
	if (0 == rank)
		printf("Time elapsed is %f seconds.\n", elapsed);
	MPI_Finalize();
	return 0;
}
