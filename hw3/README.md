# HW 3 

# Everything was run on a CIMS machine

- Program ./int_ring takes as input the number N (= 1000000 in the example below) of loops:

> mpirun -np 6 ./int_ring 1000000

The latency is estimated as elapsed/(N * size - 2) ( elapsed = time elapsed during all the communications; (N * size - 2) = total number of messages sent).
The last value, received by rank = size - 1, should be = N*np - 2.

Taking N = 1000000 and np = 6 we get ('tag' is the name of the message sent):

> The last value of tag, received from process 4 by process 5, is 5999998.
> Time elapsed is 8.093902 seconds.
> Estimated latency is 0.000001.

Running it with np = 10 on machines crunchy1 and crunchy3:

> mpirun -np 10 -host crunchy1,crunchy3 -perhost 1 ./int_ring 1000000

we get 

> 

- Program ./int_ring_largearray takes as input the number N (= 1000000 in the example below) of loops as well:

> mpirun -np 6 ./int_ring_largearray 1000000
