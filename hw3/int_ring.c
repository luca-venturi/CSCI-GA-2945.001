#include <stdio.h>
#include <mpi.h>

int main (int argc, char *argv[]) 
{	
	int size, rank, N, message = 0, i, k, tag = 0;
	MPI_Status status;

	printf( "Enter parameter N:");
   	scanf("%d", &N);

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	for (i = 0; i < N; i++) {
		if (i == 0 && rank == 0) {		
			MPI_Send(&tag, 1, MPI_INT, rank + 1, tag, MPI_COMM_WORLD);
		} else if (i > 0 && rank == 0) {
			MPI_Recv(&tag, 1, MPI_INT, size - 1, tag, MPI_COMM_WORLD, &status);
			tag += 1;
			MPI_Send(&tag, 1, MPI_INT, rank + 1, tag, MPI_COMM_WORLD);		
		}
		for (k = 1; k < size - 1; k++) {	
			if (rank == k)	{
				MPI_Recv(&tag, 1, MPI_INT, rank - 1, tag, MPI_COMM_WORLD, &status);
				tag += 1;
				MPI_Send(&tag, 1, MPI_INT, rank + 1, tag, MPI_COMM_WORLD);
			}
		}
		if (i < N -1 && rank == size - 1) {
			MPI_Recv(&tag, 1, MPI_INT, size - 2, tag, MPI_COMM_WORLD, &status);
			tag += 1;
			MPI_Send(&tag, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
		}
		else if (i = N -1 && rank == size - 1) {
			MPI_Recv(&tag, 1, MPI_INT, size - 2, tag, MPI_COMM_WORLD, &status);
			printf("The last value of tag, received by process %d of %d, is %d.\n", rank, size, tag);
		}
	}

	MPI_Finalize();
	
	return 0;
}	
