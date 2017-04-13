# homework4

The mpi_bug*.c files contain MPI bugs. The file ssort.c is a stub for parallel sample sort. 

* Jacobi Weak Scaling

 | N	| Nl | nodes | tasks | time |
 | --- | --- | --- | --- | --- |
 | 100  | 100 | 1 | 1 | 0.026222 |
 | 200 | 100 | 1 | 4 | 0.029644 |
 | 400 | 100 | 1 | 16 | 0.035653 |
 | 800 | 100 | 4 | 64 | 0.189027 |
 | 1600 | 100 | 16 | 256 | 0.156606 |
 | 3200 | 100 | 64 | 1024 | 0.395132 |
 | 6400 | 100 | 256 | 4096 | 0.459740 |

* Jacobi Strong Scaling

 | N	| Nl | nodes | tasks | time |
 | --- | --- | --- | --- | --- |
 | 6400 | 800 | 4 | 64 | 11.991474 |
 | 6400 | 400 | 16 | 256 | 2.270839 |
 | 6400 | 200 | 64 | 1024 | 1.084813 |
 | 6400 | 100 | 256 | 4096 | 0.409628 |
