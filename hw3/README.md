# HW 3 

[Everything was run on a CIMS machine]

- Program ./int_ring takes as input the number N (= 1000000 in the example below) of loops:

> mpirun -np 6 ./int_ring 1000000

The latency is estimated as elapsed/(N * size - 2) (elapsed = time elapsed during all the communications; (N * size - 2) = total number of messages sent).
The last value, received by rank = size - 1, should be = N*np - 2.

Taking N = 1000000 and np = 6 we get ('tag' is the name of the message sent):

> The last value of tag, received from process 4 by process 5, is 5999998.
> Time elapsed is 8.093902 seconds.
> Estimated latency is 0.000001.

Running it with np = 10, N =10000 on machines crunchy1 and crunchy3:

> mpirun -np 10 -host crunchy1,crunchy3 -perhost 1 ./int_ring 10000

we get 

> The last value of tag, received from process 8 by process 9, is 99998.
> Time elapsed is 100.912866 seconds.
> Estimated latency is 0.001009.

- Program ./int_ring_largearray takes as input the number N of loops as well:

> mpirun -np 6 ./int_ring_largearray 10000

The message sent is an array of 2MByte size. The bandwith is measured in Mbyte / s and estimated as (N*size-2)/(2*elapsed) (elapsed = time elapsed during all the communications; (N * size - 2) = total number of messages sent).
The last value of tag[0] ('tag' being the array message sent through), which is updated as in int_ring.c, received by rank = size - 1, should be = N*np - 2.

Taking N = 10000 and np = 6 we get:

> The last value of tag[0], received from process 4 by process 5, is 59998.
> Time elapsed is 38.029484 seconds.
> Estimated bandwitdh is of 788.835325 Mbyte/s.

