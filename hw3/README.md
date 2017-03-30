# HW 3 

- program ./int_ring takes as input the number N (=100 in the example below) of loops:

> mpirun -np 6 ./int_ring 100

The latency is estimated as elapsed/(N*size-2) (elapsed = time elapsed during all the communications; (N*size-2) = total number of messages sent)
The last value, received by rank = size - 1, should be = N*np - 2.

- Taking N = 1000000 we get:

> The last value of tag, received from process 2 by process 3, is 398.
> Time elapsed is 0.000279 seconds.
> Estimated latency is 0.000001.

