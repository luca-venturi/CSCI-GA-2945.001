/******************************************************************************
* FILE: mpi_solved2.c
* Comment: The bug was caused by the fact that an integer type message was sent by process 0 to process 1, while process 1 was expecting a float type message from process 0. The bug has been fixed by making process 1 receive an integer type message which is then converted in a float type.
******************************************************************************/
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[])
{
int numtasks, rank, tag=1, alpha, i;
float beta;
MPI_Request reqs[10];
MPI_Status stats[10];

MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

if (rank == 0) {
  if (numtasks > 2) 
    printf("Numtasks=%d. Only 2 needed. Ignoring extra...\n",numtasks);
  for (i=0; i<10; i++) {
    alpha = i*10;
    MPI_Isend(&alpha, 1, MPI_INT, 1, tag + i, MPI_COMM_WORLD, &reqs[i]);
    MPI_Wait(&reqs[i], &stats[i]);
    printf("Task %d sent = %d\n",rank,alpha);
    }
  }

if (rank == 1) {
  for (i=0; i<10; i++) {
    MPI_Irecv(&alpha, 1, MPI_INT, 0, tag + i, MPI_COMM_WORLD, &reqs[i]);
	MPI_Wait(&reqs[i], &stats[i]);	
	beta = (float) alpha;    
    printf("Task %d received = %f\n",rank,beta);
    }
  }

MPI_Finalize();

return 0;
}
