#include <stdio.h>
#include <mpi.h>
#include "util.h"

int main (int argc, char *argv[]) 
{	
	int size, rank, N, i, k;
	MPI_Status status;
	timestamp_type time1, time2;

	const int mem = 2*1024*1024;
	const int cnt = mem/sizeof(int);
	int tag[cnt];
	tag[0] = 0;
	for (i = 0; i < cnt; i++) 	
		tag[i] = i; 
			
	sscanf(argv[1], "%d", &N);

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == size - 1) 
		get_timestamp(&time1);
	/* For the following, N should be >= 2 */
	for (i = 0; i < N; i++) {
		if (i == 0 && rank == 0) {	
			MPI_Send(tag, cnt, MPI_INT, rank + 1, 199, MPI_COMM_WORLD);
		} else if (i > 0 && rank == 0) {
			MPI_Recv(tag, cnt, MPI_INT, size - 1, 199, MPI_COMM_WORLD, &status);
			tag[0] += 1;
			MPI_Send(tag, cnt, MPI_INT, rank + 1, 199, MPI_COMM_WORLD);		
		}
		for (k = 1; k < size - 1; k++) {	
			if (rank == k)	{
				MPI_Recv(tag, cnt, MPI_INT, rank - 1, 199, MPI_COMM_WORLD, &status);
				tag[0] += 1;
				MPI_Send(tag, cnt, MPI_INT, rank + 1, 199, MPI_COMM_WORLD);
			}
		}
		if (i < N -1 && rank == size - 1) {
			MPI_Recv(tag, cnt, MPI_INT, rank -1, 199, MPI_COMM_WORLD, &status);
			tag[0] += 1;			
			MPI_Send(tag, cnt, MPI_INT, 0, 199, MPI_COMM_WORLD);
		}
		else if (i == N -1 && rank == size - 1) {
			MPI_Recv(tag, cnt, MPI_INT, rank - 1, 199, MPI_COMM_WORLD, &status);
			get_timestamp(&time2);
			/* Last value should be = N*np - 2 */			
			printf("The last value of tag[0], received from process %d by process %d, is %d.\n", rank -1, rank, tag[0]);
			double elapsed = timestamp_diff_in_seconds(time1,time2);
			printf("Time elapsed is %f seconds.\n", elapsed);
			printf("Estimated bandwitdh is of %f Mbyte/s.\n", (N*size-2)/(2*elapsed));
		}
	}

	MPI_Finalize();
	
	return 0;
}	
