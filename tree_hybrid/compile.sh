#!/bin/bash
#module load openmpi/4.0.3/gnu/9.3.0 
mpicc -c sortin_n.c  -O3 -march=native -lm -fopenmp -Wall -Wextra 
mpicc tree_hybrid.c sortin_n.o -o tree_hybrid.x -O3 -lm -fopenmp -Wall -Wextra -march=native
