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
	for (i = 0; i < N; ++i)
		vec[i] = rand();
	printf("rank: %d, first entry: %d\n", rank, vec[0]);

	/* sort locally */
	qsort(vec, N, sizeof(int), compare);

	/* select local splitters, i.e., every N/P-th entry of the sorted vector */
	s = N / p;
	int *send, *recv;
	send = (int *) calloc(s, sizeof(int));
	recv = (int *) calloc(s*p, sizeof(int));
	for (i = 0; i < s; i++) 
		send[i] = vec[p*(i+1)-1]; 

	/* every processor communicates the selected entries to the root processor */
	MPI_Gather(send, s, MPI_INT, recv, s, MPI_INT, 0, MPI_COMM_WORLD);
	free(send);

	/* root processor does a sort, determinates splitters that split the data into P buckets of approximately the same size */
	int splt[p-1];
	if (0 == rank) {
		qsort(recv, s*p, sizeof(int), compare);
		for (i = 0; i < p-1; i++)
			splt[i] = recv[s*(i+1) - 1];
	}
	free(recv);

	/* root process broadcasts splitters */
	MPI_Bcast(splt, p-1, MPI_INT, 0, MPI_COMM_WORLD);

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
	free(vec);

	/* send and receive: MPI_Alltoall share with every processor how many integers it should expect, and then MPI_Send and MPI_Recv to exchange the data */
	int *bins_recv[p], bins_recv_size[p];
	MPI_Alltoall(bins_send_size, 1, MPI_INT, bins_recv_size, 1, MPI_INT, MPI_COMM_WORLD);
	for (j = 0; j < p; j++) {
		if (j != rank) {
			MPI_Send(bins_send[j], bins_send_size[j], MPI_INT, j, p*j+rank, MPI_COMM_WORLD);
		}
	}
	for (j = 0; j < p; j++) 
		bins_recv[j] = (int *) calloc(bins_recv_size[j], sizeof(int));
	for (i = 0; i < bins_recv_size[rank]; i++)		
		bins_recv[rank][i] = bins_send[rank][i];
	for (j = 0; j < p; j++) {
		if (j != rank) {
			MPI_Recv(bins_recv[j], bins_recv_size[j], MPI_INT, j, p*rank+j, MPI_COMM_WORLD, &(status[j]));
		}
	}
	int new_vec_len = 0, tmp = 0;
	for (j = 0; j < p; j++)
		new_vec_len += bins_recv_size[j];
	int * vec_new = (int *) calloc(new_vec_len, sizeof(int));
	for (i = 0; i < bins_recv_size[0]; i++) {
		vec_new[tmp] = bins_recv[0][i];
		tmp++;
	}			
	for (j = 1; j < p; j++) {
		for (i = 0; i < bins_recv_size[j]; i++) {
			vec_new[tmp] = bins_recv[j][i];
			tmp++;
		}	
	}
	for (j = 0; j < p; j++) {
		free(bins_send[j]);
		free(bins_recv[j]);
	}
	
	/* local sort */
	qsort(vec_new, new_vec_len, sizeof(int), compare);

	/* every processor writes its result to a file */
	char str[30];
	sprintf(str, "ss/p%d.txt", rank);
	FILE *f = fopen(str, "w");
	for (j = 0; j < new_vec_len; j++)
		fprintf(f, "%d\n", vec_new[j]);
	fclose(f);
	free(vec_new);
	
	MPI_Finalize();
	return 0;
}
