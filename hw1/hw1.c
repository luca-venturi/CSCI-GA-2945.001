#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"

double norm(int N, double *v);

double *jacobi(int *M, double *eps, int N, double *u, double h);

double *gauss_seidel(int *M, double *eps, int N, double *u, double h);

int main (int argc, char **argv) /* the program takes as input the number of grid points N, the method to use (0 for Jacobi, 1 for Gauss-Seidel) and the maximum number of iterations M */
{
	int i, N, *M, met;
	double *f, *u, h0, h, *eps, elapsed;
	timestamp_type time1, time2;

	if (argc != 3 && argc != 4) {
    	fprintf(stderr, "Argument not valid!\n");
    	abort();
  	}

	N = atol(argv[1]);
	met = atol(argv[2]);
	if (met != 0 && met != 1) {
    	fprintf(stderr, "Method not available!\n");
    	abort();
  	}
	h0 = N+1;
	h = (h0*h0);
	M = (int *) malloc(sizeof(int));
	if (argc == 4) {
		*M = atol(argv[3]);
	} else {
		*M = 1000; 
	}	
	eps = (double *) malloc(sizeof(double));	
	*eps = 0.0001;	
	
	u = (double *) calloc(N, sizeof(double));	
	
	get_timestamp(&time1);
	if (met == 0) {
		jacobi(M, eps, N, u, h);		
		get_timestamp(&time2);
		elapsed = timestamp_diff_in_seconds(time1,time2);
		printf("\nJacobi method for %d grid points used %d iteration to reduce the initial error by a factor of %f.\n\n", N, *M, *eps);
	} else {	
		gauss_seidel(M, eps, N, u, h);
		get_timestamp(&time2);
		elapsed = timestamp_diff_in_seconds(time1,time2);
		printf("\nGauss-Seidel method for %d grid points used %d iteration to reduce the initial error by a factor of %f.\n\n", N, *M, *eps);	
	} 
	printf("Time elapsed is %f seconds.\n\n", elapsed);

	free(M);
	free(eps);
	free(u);	

	return 0;
}

double norm(int N, double *v) 
{
	int j;
	double sum_v=0;

	for (j = 0; j < N; j++) {
		sum_v += v[j]*v[j];
	}

	return sqrt(sum_v);
}

double *jacobi(int *M, double *eps, int N, double *u, double h) 
{	
	int it = 0, max_it, j;
	double *b, *temp_u, min_eps, s, norm_b0, norm_b=1;

	min_eps = *eps;
	max_it = *M;
	norm_b0 = sqrt(N);

	temp_u = (double *) malloc(N*sizeof(double));
	b = (double *) malloc(N*sizeof(double));

	while (it < max_it && norm_b > min_eps) {
		for (j = 0; j < N; j++) {
			temp_u[j] = u[j];
		}	
		u[0] = (1+temp_u[1]*h)/(2*h);	
		for (j = 1; j < N-1; j++) {
			u[j] = (1+(temp_u[j-1]+temp_u[j+1])*h)/(2*h);
		}
		u[N-1] = (1+temp_u[N-2]*h)/(2*h);
		b[0] = -1 + h*(2*u[0]-u[1]);
		for (j = 1; j < N-1; j++) {
			b[j] = -1 + (-u[j-1]+2*u[j]-u[j+1])*h;
		}
		b[N-1] = -1 + h*(2*u[N-1]-u[N-2]);
		norm_b = norm(N, b);
		norm_b /= norm_b0;
		it++;
	}

	free(temp_u);
	free(b);
	*eps = norm_b;
	*M = it;
}

double *gauss_seidel(int *M, double *eps, int N, double *u, double h) 
{	
	int it = 0, max_it, j;
	double *b, min_eps, norm_b0, norm_b=1;

	min_eps = *eps;
	max_it = *M;
	norm_b0 = sqrt(N);

	b = (double *) malloc(N*sizeof(double));
	
	while (it < max_it && norm_b > min_eps) { 
		u[0] = (1+u[1]*h)/(2*h);	
		for (j = 1; j < N-1; j++) {
			u[j] = (1+(u[j-1]+u[j+1])*h)/(2*h);
		}
		u[N-1] = (1+u[N-2]*h)/(2*h);
		b[0] = -1 + h*(2*u[0]-u[1]);
		for (j = 1; j < N-1; j++) {
			b[j] = -1 + (-u[j-1]+2*u[j]-u[j+1])*h;
		}
		b[N-1] = -1 + h*(2*u[N-1]-u[N-2]);
		norm_b = norm(N, b);
		norm_b /= norm_b0;
		it++;
	}

	free(b);
	*eps = norm_b;
	*M = it;
}
