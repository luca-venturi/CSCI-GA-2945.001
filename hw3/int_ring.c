#include <stdio.h>
#include <mpi.h>
#include "util.h"

int main (int argc, char *argv[]) 
{	
	int size, rank, N, message = 0, i, k, tag = 0;
	MPI_Status status;
	timestamp_type time1, time2;
			
	sscanf(argv[1], "%d", &N);

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == size - 1) 
		get_timestamp(&time1);
	/* For the following, N should be >= 2 */
	for (i = 0; i < N; i++) {
		if (i == 0 && rank == 0) {	
			MPI_Send(&tag, 1, MPI_INT, rank + 1, 199, MPI_COMM_WORLD);
		} else if (i > 0 && rank == 0) {
			MPI_Recv(&tag, 1, MPI_INT, size - 1, 199, MPI_COMM_WORLD, &status);
			tag += 1;
			MPI_Send(&tag, 1, MPI_INT, rank + 1, 199, MPI_COMM_WORLD);		
		}
		for (k = 1; k < size - 1; k++) {	
			if (rank == k)	{
				MPI_Recv(&tag, 1, MPI_INT, rank - 1, 199, MPI_COMM_WORLD, &status);
				tag += 1;
				MPI_Send(&tag, 1, MPI_INT, rank + 1, 199, MPI_COMM_WORLD);
			}
		}
		if (i < N -1 && rank == size - 1) {
			MPI_Recv(&tag, 1, MPI_INT, rank -1, 199, MPI_COMM_WORLD, &status);
			tag += 1;
			MPI_Send(&tag, 1, MPI_INT, 0, 199, MPI_COMM_WORLD);
		}
		else if (i == N -1 && rank == size - 1) {
			MPI_Recv(&tag, 1, MPI_INT, rank - 1, 199, MPI_COMM_WORLD, &status);
			get_timestamp(&time2);
			/* Last value should be = N*np - 2 */			
			printf("The last value of tag, received from process %d by process %d, is %d.\n", rank -1, rank, tag);
			double elapsed = timestamp_diff_in_seconds(time1,time2);
			printf("Time elapsed is %f seconds.\n", elapsed);
			printf("Estimated latency is %f.\n", elapsed/(N*size-2));
		}
	}

	MPI_Finalize();
	
	return 0;
}	
