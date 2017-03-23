#include <stdio.h>
#include <math.h>
#include "util.h"
#include <string.h>

/* compute global residual, assuming ghost values are updated */
double compute_residual(double *u, int Ntot, int Ntotsq, double invhsq)
{
	int i;
	double tmp, res = 0.0;
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

	/* timing */
	timestamp_type time1, time2;
	get_timestamp(&time1);

	/* Allocation of vectors, including ghost points */
	double * u    = (double *) calloc(sizeof(double), Ntotsq);
	double * unew = (double *) calloc(sizeof(double), Ntotsq);	
	double h = 1.0 / (N + 1); 
	double hsq = h * h;
	double invhsq = 1./hsq;
	double res, res0, tol = 1e-5;

	/* initial residual */
	res0 = N;
	res = res0;

	for (iter = 0; iter < max_iters && res/res0 > tol; iter++) {

    	/* Jacobi step for all the inner points */
    	for (i = Ntot+1; i <= Ntotsq-Ntot-1; i++) {
			if ((i % Ntot) != 0 && (i % Ntot) != Ntot-1) {
				unew[i] = 0.25 * (hsq + u[i-1] + u[i+1] + u[i+Ntot] + u[i-Ntot]);			
			}	
		}

		/* copy new_u onto u */
		double *utemp;
		utemp = u;	
		u = unew;
		unew = utemp;
		
		if (0 == (iter % 10)) {
      		res = compute_residual(u, Ntot, Ntotsq, invhsq);  	
		}
	}

	printf("\nIter: %d. Residual: %g\n", iter, res/res0);

	/* Clean up */
	free(u);
	free(unew);

	/* timing */
	get_timestamp(&time2);
	double elapsed = timestamp_diff_in_seconds(time1,time2);
	printf("\nTime elapsed is %f seconds.\n\n", elapsed);

	return 0;
}
