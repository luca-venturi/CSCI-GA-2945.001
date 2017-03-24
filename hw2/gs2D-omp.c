#include <stdio.h>
#include <math.h>
#include "util.h"
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* compute global residual for black labels, assuming values are updated */
double compute_residual_b(double *u_b, double *u_r, int Ntot, int Nb_s, int Nb_e, int Nb_r1, double Nb_r2, int Nb_r3, int Nb_flag, int Nb_a, int Nb_b, double invhsq)
{
	int i, i_mod;
	double tmp, res = 0.0;
#pragma omp parallel for default(none) shared(u_b, u_r, Ntot, Nb_s, Nb_e, Nb_r1, Nb_r2, Nb_r3, Nb_flag, Nb_a, Nb_b, invhsq) private(i, i_mod, tmp) reduction(+:res)
	for (i = Nb_s; i <= Nb_e; i++) {
		i_mod = (i % Ntot); 
		if (i_mod > Nb_r1 && i_mod < Nb_r2) {
			tmp = (1.0 + (u_r[i-1] + u_r[i] + u_r[i+Nb_a] + u_r[i-Nb_b] - 4.0*u_b[i]) * invhsq);
			res += tmp * tmp;
		}
		else if (i_mod > Nb_r2 && i_mod < Nb_r3) {
			tmp = (1.0 + (u_r[i+Nb_flag] + u_r[i] + u_r[i+Nb_a] + u_r[i-Nb_b] - 4.0*u_b[i]) * invhsq);
			res += tmp * tmp;
		}	
	}
	return res;
}

/* compute global residual for red labels, assuming values are updated */
double compute_residual_r(double *u_r, double *u_b, int Ntot, int Nr_s, int Nr_e, int Nr_r1, int Nr_r2, int Nr_r3, int Nr_flag, int Nr_a, int Nr_b, double invhsq)
{
	int i, i_mod;
	double tmp, res = 0.0;
#pragma omp parallel for default(none) shared(u_r, u_b, Ntot, Nr_s, Nr_e, Nr_r1, Nr_r2, Nr_r3, Nr_flag, Nr_a, Nr_b, invhsq) private(i, i_mod, tmp) reduction(+:res)
	for (i = Nr_s; i <= Nr_e; i++) {
		i_mod = (i % Ntot);			
		if (i_mod < Nr_r1) {
			tmp = (1.0 + (u_b[i] + u_b[i+1] + u_b[i+Nr_a] + u_b[i-Nr_b] - 4.0*u_r[i]) * invhsq);
			res += tmp * tmp;
		}	
		else if (i_mod > Nr_r2 && i_mod < Nr_r3) {
			tmp = (1.0 + (u_b[i+Nr_flag] + u_b[i] + u_b[i+Nr_a] + u_b[i-Nr_b] - 4.0*u_r[i]) * invhsq);
			res += tmp * tmp;
		}
	}
	return res;
}

int main(int argc, char * argv[])
{
	int i, i_mod, N, iter, max_iters;

	sscanf(argv[1], "%d", &N);
	int Ntot = N + 2;
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
		
    /* number of black and red variables */
	int Nb = N*(Ntot+2)*0.5 + 2.5;
	int Nr = Ntot*Ntot*0.5;

	/* timing */
	timestamp_type time1, time2;
	get_timestamp(&time1);

	/* Allocation of vectors, including ghost points */
	double * u_b    = (double *) calloc(sizeof(double), Nb);
	double * u_r    = (double *) calloc(sizeof(double), Nr);
	double h = 1.0 / (N + 1); 
	double hsq = h * h;
	double invhsq = 1./hsq;
	double res, res0, tol = 1e-5, tmp;

	/* initial residual */
	res0 = N;
	res = res0;

	/* additional labeling variables */
	int Nb_s = (Ntot+1)*0.5;
	int Nb_e = (Ntot)*(Ntot-1)*0.5 - 1;
	int Nb_r1 = 0;
	double Nb_r2 = (N+1)*0.5;
	int Nb_r3 = 2*Nb_r2 + N%2;
	int Nb_flag = 1 - 2*(N%2);
	int Nb_a = Ntot*0.5;
	int Nb_b = Nb_a + N%2;

	int Nr_s = Ntot*0.5 + 1;
	int Nr_e = (Ntot)*(Ntot-1)*0.5 - 2 - (N == 1);
	int Nr_r1 = (N + (N%2))*0.5;
	int Nr_r2 = Nr_r1 + 1 - N%2;
	int Nr_r3 = 2*Nr_r2;
	int Nr_flag = -1 + 2*(N%2);
	int Nr_a = (Ntot+1)*0.5;
	int Nr_b = Nr_a - (N%2);

	for (iter = 0; iter < max_iters && res/res0 > tol; iter++) {

    	/* GS step for all the black points */
#pragma omp parallel for default(none) shared(Ntot, Nb_s, Nb_e, Nb_r1, Nb_r2, Nb_r3, Nb_flag, Nb_a, Nb_b, u_b, u_r, hsq) private(i, i_mod) schedule(dynamic,10)
    	for (i = Nb_s; i <= Nb_e; i++) {
			i_mod = (i % Ntot); 
			if (i_mod > Nb_r1 && i_mod < Nb_r2) {
				u_b[i] = 0.25 * (hsq + u_r[i-1] + u_r[i] + u_r[i+Nb_a] + u_r[i-Nb_b]); 
			}
			else if (i_mod > Nb_r2 && i_mod < Nb_r3) {
				u_b[i] = 0.25 * (hsq + u_r[i+Nb_flag] + u_r[i] + u_r[i+Nb_a] + u_r[i-Nb_b]); 
			}	
		}

		/* GS step for all the red points */
#pragma omp parallel for default(none) shared(Ntot, Nr_s, Nr_e, Nr_r1, Nr_r2, Nr_r3, Nr_flag, Nr_a, Nr_b, u_b, u_r, hsq) private(i, i_mod) schedule(dynamic,10)
    	for (i = Nr_s; i <= Nr_e; i++) {
			i_mod = (i % Ntot);			
			if (i_mod < Nr_r1) {
				u_r[i] = 0.25 * (hsq + u_b[i+1] + u_b[i] + u_b[i+Nr_a] + u_b[i-Nr_b]);
			}	
			else if (i_mod > Nr_r2 && i_mod < Nr_r3) {
				u_r[i] = 0.25 * (hsq + u_b[i+Nr_flag] + u_b[i] + u_b[i+Nr_a] + u_b[i-Nr_b]);
			}
		}
		
		if (0 == (iter % 10)) {
			tmp = compute_residual_b(u_b, u_r, Ntot, Nb_s, Nb_e, Nb_r1, Nb_r2, Nb_r3, Nb_flag, Nb_a, Nb_b, invhsq);
			tmp += compute_residual_r(u_r, u_b, Ntot, Nr_s, Nr_e, Nr_r1, Nr_r2, Nr_r3, Nr_flag, Nr_a, Nr_b, invhsq);
			res = sqrt(tmp);
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
