/* Parallel sample sort */
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>

static int compare(const void *a, const void *b)
{
	int *da = (int *)a;
	int *db = (int *)b;

	if (*da > *db)
		return 1;
	else if (*da < *db)
		return -1;
	else
		return 0;
}

int main( int argc, char *argv[])
{
	int rank, p;
	int i, j, N, s;
	int *vec;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Status status[p-1]; 

	/* Number of random numbers per processor */
	sscanf(argv[1], "%d", &N);

	vec = (int *) calloc(N, sizeof(int));
	/* seed random number generator differently on every core */
	srand((unsigned int) (rank + 393919));

	/* fill vector with random integers */
	for (i = 0; i < N; ++i) {
		vec[i] = rand();
	}
	printf("rank: %d, first entry: %d\n", rank, vec[0]);

	/* sort locally */
	qsort(vec, N, sizeof(int), compare);

    /* randomly sample s entries from vector or select local splitters, i.e., every N/P-th entry of the sorted vector */
	s = (int)(N / p);
	int *send, *recv;
	send = (int *) calloc(s, sizeof(int));
	if (0 == rank) {
		recv = (int *) calloc(s*p, sizeof(int));
	}
	for (i = 0; i < s; i++) 
		send[i] = vec[p*(i+1)-1]; 

	/* every processor communicates the selected entries to the root processor */
	MPI_Gather(&(send[0]), s, MPI_INT, &(recv[0]), s, MPI_INT, 0, MPI_COMM_WORLD);

	/* root processor does a sort, determinates splitters that split the data into P buckets of approximately the same size */
	int splt[p-1];
	if (0 == rank) {
		qsort(recv, s*p, compare);
		for (i = 0; i < p-1; i++)
			splt[i] = recv[s*(i+1) - 1];
	}

	/* root process broadcasts splitters */
	MPI_Scatter(&(splt[0]), p-1, MPI_INT, 0, MPI_COMM_WORLD);

	/* every processor uses the obtained splitters to decide which integers need to be sent to which other processor (local bins) */
	int *bins_send[p], bins_send_size[p], loc_splt[p-1];
	i = 0;
	while (vec[i] <= splt[0])
		i++;
	loc_splt[0] = i;
	for (j = 1; j < p-1; j++) {
		while (splt[j-1] < vec[i] && vec[i] <= splt[j])
			i++;
		loc_splt[j] = i;
	}
	bins_send_size[0] = loc_splt[0];
	for (j = 1; j < p-1; j++)
		bins_send_size[j] = loc_splt[j] - loc_splt[j-1];
	bins_send_size[p-1] = N - loc_splt[p-2];
	bins_send[0] = (int *) calloc(bins_send_size[0], sizeof(int));
	for (i = 0; i < bins_send_size[0]; i++)
		bins_send[0][i] = vec[i];
	for (j = 1; j < p; j++) {
		bins_send[j] = (int *) calloc(bins_send_size[j], sizeof(int));
		for (i = 0; i < bins_send_size[j]; i++)
			bins_send[j][i] = vec[loc_splt[j-1]+i];
	}

	/* send and receive: either you use MPI_AlltoallV, or (and that might be easier), use an MPI_Alltoall to share with every processor how many integers it should expect, and then use MPI_Send and MPI_Recv to exchange the data */
	int *bins_recv[p], bins_recv_size[p];
	MPI_Alltoall(bins_send_size, 1, MPI_INT, bins_recv_size, 1, MPI_INT, MPI_COMM_WORLD);
	for (j = 0; j < p; j++) {
		if (j != rank) {
			MPI_Send(bins_send[j], bins_send_size[j], MPI_INT, j, p*j+rank, MPI_COMM_WORLD);
		}
	}
	for (j = 0; j < p; j++) 
		bins_recv[j] = (int *) calloc(bins_recv_size[j], sizeof(int));
	bins_recv[rank] = bins_send[rank];
	for (j = 0; j < p; j++) {
		if (j != rank) {
			MPI_Recv(bins_recv[j], bins_recv_size[j], MPI_INT, p*rank+j, tag, MPI_COMM_WORLD, &(status[j]));
		}
	}
	for (i = 0; i < bins_recv_size[0]; i++)
		vec[i] = bins_recv[0][i];
	for (j = 1; j < p; j++) {
		for (i = 0; i < bins_recv_size[j]; i++)
			vec[bins_recv_size[j-1]+i] = bins_recv[j][i];
	}
	
	
	/* do a local sort */
	qsort(vec, N, sizeof(int), compare);

	/* every processor writes its result to a file */

	free(vec);
	MPI_Finalize();
	return 0;
}
