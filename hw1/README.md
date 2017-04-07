## HW1

### Jacobi/Gauss-Seidel for 1D diffusion equation

Program ./hw1 takes as input the number N (= 100 in the example below) of discretization points and the maximum number of iterations M (= 1000 in the example below):

> ./hw1 100 1000

##### Timings

Running on a CIMS machine I got:

* Jacobi

 | N	| -O0 | -O3	|
 | --- | --- | --- |
 | 1000  | 0.0311 | 0.0111 |
 | 100000 | 1.8277 | 0.6680 |

* Gauss-Seidel

 | N	| -O0 | -O3	|
 | --- | --- | --- |
 | 1000  | 0.0381 | 0.0210 |
 | 100000 | 2.1496 | 1.2982 |

---

Running on my machine I got:

* Jacobi

 | N	| -O0 | -O3	|
 | --- | --- | --- |
 | 1000  | 0.0408 | 0.0151 |
 | 100000 | 2.5506 | 0.8218 |

* Gauss-Seidel

 | N	| -O0 | -O3	|
 | --- | --- | --- |
 | 1000  | 0.0344 | 0.0354 |
 | 100000 | 3.4644 | 2.1768 |

My machine (PC) has a processor Intel-Core i3-4010U CPU @ 1.70 GHz with 2 cores and a 4
GB total RAM. OS: is Linux Ubuntu 16.04 (64-bit).
