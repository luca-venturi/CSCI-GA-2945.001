#include <stdio.h>
#include <math.h>
#include "util.h"
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* compute global residual, assuming ghost values are updated */
double compute_residual(double *u, int Ntot, int Ntotsq, double invhsq)
{
	int i;
	double tmp, res = 0.0;
#pragma omp parallel for default(none) shared(u, Ntot, Ntotsq, invhsq) private(i, tmp) reduction(+:res)
	for (i = Ntot+1; i <= Ntotsq-Ntot-1; i++) {
		if ((i % Ntot) != 0 && (i % Ntot) != Ntot-1) {
			tmp = (1.0 + (u[i-1] + u[i+1] + u[i+Ntot] + u[i-Ntot] - 4.0*u[i]) * invhsq);
			res += tmp * tmp;
		}
	}
	return sqrt(res);
}

int main(int argc, char * argv[])
{
	int i, N, iter, max_iters;

	sscanf(argv[1], "%d", &N);
	int Ntot = N + 2;
	int Ntotsq = Ntot * Ntot;
	sscanf(argv[2], "%d", &max_iters);
	#pragma omp parallel
	{
#ifdef _OPENMP
		int my_threadnum = omp_get_thread_num();
		int numthreads = omp_get_num_threads();
#else
		int my_threadnum = 0;
		int numthreads = 1;
#endif
		printf("Hello, I'm thread %d out of %d\n", my_threadnum, numthreads);
	}
		
	int Nb = N*(Ntot+2)*0.5 + 2.5;
	int Nr = Ntot*Ntot*0.5;

	/* timing */
	timestamp_type time1, time2;
	get_timestamp(&time1);

	/* Allocation of vectors, including left and right ghost points */
	double * u_b    = (double *) calloc(sizeof(double), Nb);
	double * u_r    = (double *) calloc(sizeof(double), Nr);
	double h = 1.0 / (N + 1); 
	double hsq = h * h;
	double invhsq = 1./hsq;
	double res, res0, tol = 1e-5;

	/* initial residual */
	res0 = N;
	res = res0;

	int Nb_s = (Ntot+1)*0.5;
	int Nb_e = (Ntot)*(Ntot-1)*0.5 - 1;
	int Nb_mod = Ntot;
	int Nb_r1 = 0;
	int Nb_r2 = (N+1)*0.5*(2-N%2);
	int Nb_a = Ntot*0.5;
	int Nb_b = Nb_a + N%2;
	int Nb_l;

	int Nr_s = Ntot*0.5 + 1;
	int Nr_e = (Ntot)*(Ntot-1)*0.5 - 2 - (N == 1);
	int Nr_mod = Ntot;
	int Nr_r1 = (N + (N%2))*0.5;
	int Nr_r2 = N*(N%2+1)*0.5 + 1;
	int Nr_a = 0;
	int Nr_b = Ntot*0.5;
	int Nr_l;
	printf("%d", Nr_b);

	for (iter = 0; iter < max_iters && res/res0 > tol; iter++) {

    	/* Jacobi step for all the black points */
    	for (i = Nb_s; i <= Nb_e; i++) {
			if ((i % Nb_mod) != Nb_r1 && (i % Nb_mod) != Nb_r2) {
				u_b[i] = 0.25 * (hsq + u_r[i+Nb_l] + u_r[i] + u_r[i+Nb_a] + u_r[i-Nb_b]); // left !!!
			}	
		}

		/* Jacobi step for all the red points */
    	for (i = Nr_s; i <= Nr_e; i++) {
			if ((i % Nr_mod) != Nr_r1 && (i % Nr_mod) != Nr_r2) {
				u_r[i] = 0.25 * (hsq + u_b[i+Nr_l] + u_b[i] + u_b[i+Nr_a] + u_b[i-Nr_b]); // above - right - left !!!			
			}	
		}
		
		if (0 == (iter % 10)) {
      		res = compute_residual(u_r, Ntot, Ntotsq, invhsq);  	////////////////////////////////////////////////////////
		}
	}

	printf("\nIter: %d. Residual: %g\n", iter, res/res0);

	/* Clean up */
	free(u_b);
	free(u_r);

	/* timing */
	get_timestamp(&time2);
	double elapsed = timestamp_diff_in_seconds(time1,time2);
	printf("\nTime elapsed is %f seconds.\n\n", elapsed);
	
	return 0;
}
