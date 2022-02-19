gcc -c sortin_n.c  -O3 -march=native -lm -fopenmp -Wall -Wextra 
gcc tree_sorting_network.c sortin_n.o -o tree_sn_32.x -O3 -lm -fopenmp -Wall -Wextra -march=native
