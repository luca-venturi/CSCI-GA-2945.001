#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"

double norm(int N, double *v);

void *jacobi(int *M, double *eps, int N, double *u, double h);

int main (int argc, char **argv) /* the program takes as input the number of grid points N and the maximum number of iterations M */
{
	int N, *M;
	double *u, h, *eps, elapsed;
	timestamp_type time1, time2;

	if (argc != 3 && argc != 2) {
    	fprintf(stderr, "Argument not valid!\n");
    	abort();
  	}

	N = atol(argv[1]);
	h = 1/(double)((N+1)*(N+1));
	M = (int *) malloc(sizeof(int));
	if (argc == 3) {
		*M = atol(argv[2]);
	} else {
		*M = 1000; 
	}	
	eps = (double *) malloc(sizeof(double));	
	*eps = 0.0001;	
	
	u = (double *) calloc((N+2)*(N+2), sizeof(double));	
	
	get_timestamp(&time1);
	jacobi(M, eps, N, u, h);		
	get_timestamp(&time2);
	elapsed = timestamp_diff_in_seconds(time1,time2);
	printf("\nJacobi method for %d grid points used %d iteration to reduce the initial error by a factor of %f.\n\n", N, *M, *eps);
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

void *jacobi(int *M, double *eps, int N, double *u, double h) 
{	
	int it = 0, max_it, i, j, N2;
	double *b, *temp_u, min_eps, norm_b0, norm_b=1;

	min_eps = *eps;
	max_it = *M;
	norm_b0 = N;
	N2 = (N+2)*(N+2);

	temp_u = (double *) calloc(N2, sizeof(double));
	b = (double *) calloc(N2, sizeof(double));

	while (it < max_it && norm_b > min_eps) {
		for (i = 0; i<N2; i++) 
			temp_u[i] = u[i];
		
		for (i = 1; i < N+1; i++) {
			for (j = 1; j < N+1; j++) 
				u[N*i+j] = (h - temp_u[N*(i-1)+j] - temp_u[N*(i+1)+j] - temp_u[N*i+j+1] - temp_u[N*i+j-1])/4.0;
			/* for (j = 1; j < N+1; j++) {
				u[N*i+j] = (h - temp_u[N*(i-1)+j] - temp_u[N*(i+1)+j] - temp_u[N*i+j+1] - temp_u[N*i+j-1])/4.0;
				printf("%f", u[N*i+j]); } printf("\n"); */
		}
		for (i = 1; i < N+1; i++) {
			for (j = 1; j < N+1; j++)
				b[N*i+j] = 1.0 + (u[N*(i-1)+j] + u[N*(i+1)+j] + u[N*i+j+1] + u[N*i+j-1] - 4.0*u[N*i+j])/h;
		}
			
		norm_b = norm(N2, b);
		norm_b /= norm_b0;
		printf("%f\n", norm_b); /**/
		it++;
	}

	free(temp_u);
	free(b);
	*eps = norm_b;
	*M = it;
}
