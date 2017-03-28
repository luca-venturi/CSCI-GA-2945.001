#include <stdio.h>
#include <mpi.h>

int main (int argc, char *argv[]) 
{	
	int size, rank, N, message = 0, i, k;

	printf( "Enter parameter N:");
   	scanf("%d", &N);

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	for (i = 0; i < N; i++) {
		for (k = 0; k < size; k++) {	
			if (rank == k)	{
				
			}
		}
	}

	MPI_Finalize();
	
	return 0;
}	
