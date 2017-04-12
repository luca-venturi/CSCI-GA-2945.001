/* MPI-parallel Jacobi smoothing to solve - (u_xx + u_yy) = f
 * Global vector has (N+2)*(N+2) unknowns, each processor (p total) works with its
 * part, which has (lN+2)*(lN+2) unknowns, where lN = N / sqrt(p), p = 4**j.
 * Author: Luca Venturi
 */
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "util.h"
#include <string.h>

/* compute global residual, assuming ghost values are updated */
double compute_residual(double *lu, int lN, int lNtot, double invhsq)
{
	int i, j, tmp1;
	double tmp, gres = 0.0, lres = 0.0;

	for (i = 1; i <= lN; i++) {
		for (j = 1; j <= lN; j++) {
			tmp1 = i + lNtot*j;
			tmp = (1.0 + (lu[tmp1-1] + lu[tmp1+1] + lu[tmp1+lNtot] + lu[tmp1-lNtot] - 4.0*lu[tmp1])*invhsq);
			lres += tmp * tmp;
		}	
	}
	/* use allreduce for convenience; a reduce would also be sufficient */
	MPI_Allreduce(&lres, &gres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return sqrt(gres);
}


int main(int argc, char * argv[])
{
	int mpirank, i, p, lp, tmp, j, N, lN, lNtot, lNtotsq, iter, max_iters;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	/* get name of host running MPI process */
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(processor_name, &name_len);
	printf("Rank %d/%d running on %s.\n", mpirank, p, processor_name);

	sscanf(argv[1], "%d", &N);
	sscanf(argv[2], "%d", &max_iters);

	/* get j = log_4 (p) */
	j = 0;
	tmp = p;
	while (tmp > 1) {
		tmp /= 4;
		j++;
	}
	if (tmp != 1) {
		printf("Exiting. p must be a power of 4\n");
		MPI_Abort(MPI_COMM_WORLD, 0);
	}
	lp = 1;
	for (i = 1; i <= j; i++)
		lp *= 2;
	/* compute number of unknowns handled by each process */
	lN = N / lp;
	lNtot = lN + 2;
	lNtotsq = lNtot * lNtot; 
	if ((N % lp != 0) && mpirank == 0 ) {
		printf("N: %d, local N: %d\n", N, lN);
		printf("Exiting. N must be a multiple of sqrt(p)\n");
		MPI_Abort(MPI_COMM_WORLD, 0);
	}
	MPI_Status statusa, statusb, statusl, statusr;
	/* timing */
	MPI_Barrier(MPI_COMM_WORLD);
	timestamp_type time1, time2;
	get_timestamp(&time1);

	/* Allocation of vectors, including left/upper and right/lower ghost points */
	double * lu    = (double *) calloc(sizeof(double), lNtotsq);
	double * lunew = (double *) calloc(sizeof(double), lNtotsq);
	double * lutemp;
	double *lughost_left, *lughost_right;
	if (mpirank % lp != lp -1) {
		lughost_left = (double *) malloc(sizeof(double)*lN);
	}
	if (mpirank % lp != 0) {
		lughost_right = (double *) malloc(sizeof(double)*lN);
	}

	double h = 1.0 / (N + 1);
	double hsq = h * h;
	double invhsq = 1./hsq;
	double gres, gres0, tol = 1e-5;

	/* initial residual */
	gres0 = compute_residual(lu, lN, lNtot, invhsq); /**/
	gres = gres0; /**/

	for (iter = 0; iter < max_iters && gres/gres0 > tol; iter++) {

		/* Jacobi step for local points */
		for (i = 1; i <= lN; i++) {
			for (j = 1; j <= lN; j++) {
				tmp = i + lNtot*j;
				lunew[tmp]  = 0.25 * (hsq + lu[tmp - 1] + lu[tmp + 1] + lu[tmp + lNtot] + lu[tmp - lNtot]);
			}
		}

		/* communicate ghost values above - below */
		/* send: above -> below; receive: above <- below */
		if (mpirank >= lp) {
			MPI_Send(&(lunew[lNtot+1]), lN, MPI_DOUBLE, mpirank-lp, 124, MPI_COMM_WORLD);
			MPI_Recv(&(lunew[1]), lN, MPI_DOUBLE, mpirank-lp, 123, MPI_COMM_WORLD, &statusa);
		}
		/* send: below -> above; receive: below <- above */
		if (mpirank < p - lp) {
			MPI_Send(&(lunew[lNtotsq-lNtot-lN-1]), lN, MPI_DOUBLE, mpirank+lp, 123, MPI_COMM_WORLD);
			MPI_Recv(&(lunew[lNtotsq-lN-1]), lN, MPI_DOUBLE, mpirank+lp, 124, MPI_COMM_WORLD, &statusb);
		}

		/* communicate ghost values right - left */
		/* send: left -> right; receive: left <- right */
		if (mpirank % lp != lp -1) {
			for (i = 0; i < lN; i++)
				lughost_left[i] = lunew[lNtot*(i+1)+lN];
			MPI_Send(&(lughost_left[0]), lN, MPI_DOUBLE, mpirank+1, 128, MPI_COMM_WORLD);
			MPI_Recv(&(lughost_left[0]), lN, MPI_DOUBLE, mpirank+1, 129, MPI_COMM_WORLD, &statusl);
			for (i = 0; i < lN; i++)
				lunew[lNtot*(i+1)+lN+1] = lughost_left[i];
		}
		/* send: right -> left; receive: right <- left */
		if (mpirank % lp != 0) {
			for (i = 0; i < lN; i++)
				lughost_right[i] = lunew[lNtot*(i+1)+1];
			MPI_Send(&(lughost_right[0]), lN, MPI_DOUBLE, mpirank-1, 129, MPI_COMM_WORLD);
			MPI_Recv(&(lughost_right[0]), lN, MPI_DOUBLE, mpirank-1, 128, MPI_COMM_WORLD, &statusr);
			for (i = 0; i < lN; i++)
				lunew[lNtot*(i+1)] = lughost_right[i];
		}


		/* copy newu to u using pointer flipping */
		lutemp = lu; lu = lunew; lunew = lutemp;
		if (0 == (iter % 10)) {
			gres = compute_residual(lu, lN, lNtot, invhsq); 
			if (0 == mpirank) {
				printf("Iter %d: Residual: %g\n", iter, gres/gres0);
			}
		}
	}

	/* Clean up */
	free(lu);
	free(lunew);
	if (mpirank % lp != lp -1) {
		free(lughost_left);
	}
	if (mpirank % lp != 0) {
		free(lughost_right);
	}

	/* timing */
	MPI_Barrier(MPI_COMM_WORLD);
	get_timestamp(&time2);
	double elapsed = timestamp_diff_in_seconds(time1,time2);
	if (0 == mpirank) {
		printf("Time elapsed is %f seconds.\n", elapsed);
	}
	MPI_Finalize();
	return 0;
}
