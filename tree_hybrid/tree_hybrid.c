#define _GNU_SOURCE // needed to get binding information
#include <sched.h>

#include "sortin_n.h"
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define MASTER 0

knode *global_node_address = NULL;
int global_slave = 0;

void init_point(kpoint *dati, unsigned int n) {
    //#pragma omp parallel for
    for (unsigned int i = 0; i < n; ++i) {
        dati[i].coord[0] = drand48();
        dati[i].coord[1] = drand48();
    }
}

void swap(kpoint *restrict a, kpoint *restrict b) {
    kpoint t = *a;
    a->coord[0] = b->coord[0];
    a->coord[1] = b->coord[1];

    b->coord[0] = t.coord[0];
    b->coord[1] = t.coord[1];
}

void view_tree(knode *dati) {
    if ((dati->right == NULL) && (dati->left == NULL)) {
        return;
    }

    if (dati->right == NULL) {
        kpoint point = dati->split;
        knode *left = dati->left;
        kpoint l_point = left->split;
        printf(" \"%f,%f\" -- \"%f,%f\" [label=%c] \n", point.coord[0], point.coord[1], l_point.coord[0], l_point.coord[1], dati->axis + 88);

        view_tree(left);
        return;
    }

    if (dati->left == NULL) {
        kpoint point = dati->split;
        knode *right = dati->right;
        kpoint r_point = right->split;
        printf(" \"%f,%f\" -- \"%f,%f\" [label=%c]\n", point.coord[0], point.coord[1], r_point.coord[0], r_point.coord[1], dati->axis + 88);

        view_tree(right);
        return;
    }

    kpoint point = dati->split;
    knode *right = dati->right;
    knode *left = dati->left;

    kpoint r_point = right->split;
    kpoint l_point = left->split;

    printf(" \"%f,%f\" -- \"%f,%f\" [label=%c] \n", point.coord[0], point.coord[1], r_point.coord[0], r_point.coord[1], dati->axis + 88);
    printf(" \"%f,%f\" -- \"%f,%f\" [label=%c] \n", point.coord[0], point.coord[1], l_point.coord[0], l_point.coord[1], dati->axis + 88);

    view_tree(left);
    view_tree(right);
    return;
}

unsigned int three_way_partition(kpoint *dati, TYPE midvalue, unsigned int npoint, unsigned short int axis) {
    unsigned int i = 0;
    unsigned int j = 0;
    unsigned int k = npoint - 1;
    while (j <= k) {
        if (dati[j].coord[axis] < midvalue) {
            swap(dati + i, dati + j);
            ++i;
            ++j;
        } else if (dati[j].coord[axis] > midvalue) {
            swap(dati + j, dati + k);
            --k;
        } else {
            ++j;
        }
    }
    // printf("j %d, k %d\n",j,k);
    return j;
}

TYPE find_mean(kpoint *dati, unsigned int npoint, unsigned short int axis) {
    TYPE min = dati[0].coord[axis];
    TYPE max = dati[0].coord[axis];
    for (unsigned int i = 0; i < npoint; ++i) {
        if (dati[i].coord[axis] > max) {
            max = dati[i].coord[axis];
            continue;
        }
        if (dati[i].coord[axis] < min) {
            min = dati[i].coord[axis];
        }
    }

    return (min + max) / 2;
}

int find_pivot(kpoint *dati, TYPE mean, TYPE npoint, int axis) {
    for (int i = 0; i < npoint; ++i) {
        if (dati[i].coord[axis] == mean)
            return (i);
    }
    return -1;
}

TYPE find_kth(kpoint *dati, unsigned int k, unsigned short int axis, unsigned int npoints) {

    TYPE pivot = dati[0].coord[axis];

    unsigned int split_point = three_way_partition(dati, pivot, npoints, axis);

    kpoint *l = dati;
    kpoint *r = dati + split_point;

    unsigned int size_r = split_point - 1;
    unsigned int size_l = npoints - split_point;

    if (k == size_l + 1)
        return pivot;
    if (k < size_l + 1)
        return find_kth(r, k, axis, size_l);
    if (k > size_l + 1)
        return find_kth(l, k - (size_l + 1), axis, size_r);
    return -1;
}

knode *buildtree(kpoint *dati, int ndim, int axis, unsigned int size, knode *base_address, knode *memory_offset) {

    int myaxis = (axis + 1) % ndim;
    knode *this_node = NULL;
#pragma omp atomic capture
    this_node = global_node_address++;

    this_node->axis = myaxis;

    if (size == 1) { // only one node, null leaf and stop here
        this_node->left = NULL;
        this_node->right = NULL;
        this_node->split.coord[0] = dati->coord[0];
        this_node->split.coord[1] = dati->coord[1];
        return this_node;
    }

    if (size == 2) { // onle 2 node, one null leaf and build another one node

        this_node->split.coord[0] = dati->coord[0];
        this_node->split.coord[1] = dati->coord[1];
        knode *res = buildtree(++dati, 2, myaxis, --size, base_address, memory_offset);
        this_node->right = (res - base_address) + memory_offset;
        this_node->left = NULL;
        /* this_node->right = buildtree(++dati, 2, myaxis, --size,base_address,memory_offset); */
        /* knode* test=(this_node->right-base_address)+memory_offset; */
        /* this_node->right=test; */
#ifdef DEBUG
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        printf("%p right [%d]\n", this_node->right, rank);
#endif
        return this_node;
    }

    unsigned int median = ((size % 2) == 0) ? size / 2 : size / 2 + 1;

    if (size > 15) {
        TYPE median_value = find_kth(dati, median, myaxis, size);
        // median = three_way_partition(dati, median_value, size, myaxis);
    } else {
        (*sortin_network[size - 1])(dati, myaxis);
    }
    --median;

    this_node->split.coord[0] = dati[median].coord[0];
    this_node->split.coord[1] = dati[median].coord[1];

    kpoint *l_slice = dati;
    kpoint *r_slice = dati + median + 1;
    unsigned int l_size = median;
    unsigned int r_size = size - median - 1;

#pragma omp task shared(dati) firstprivate(myaxis, median, size)
    {
        knode *res = (buildtree(l_slice, 2, myaxis, l_size, base_address, memory_offset));
        this_node->left = (res - base_address) + memory_offset;

#ifdef DEBUG
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        printf("%p left [%d]\n", this_node->left, rank);
#endif
    }
#pragma omp task shared(dati) firstprivate(myaxis, median, size)
    {
        knode *res = (buildtree(r_slice, 2, myaxis, r_size, base_address, memory_offset));
        this_node->right = (res - base_address) + memory_offset;


#ifdef DEBUG
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        printf("%p right [%d]\n", this_node->right, rank);
#endif
    }
    return this_node;
}

knode *buildtree_master(kpoint *dati, int ndim, int axis, unsigned int size, int dept, int max_dept) {

    int myaxis = (axis + 1) % ndim;
    knode *this_node = malloc(sizeof(knode));
    /* #pragma omp atomic capture */
    /*   this_node=global_node_address++; */
    this_node->axis = myaxis;
    if (size == 1) { // only one node, null leaf and stop here
        this_node->left = NULL;
        this_node->right = NULL;
        this_node->split.coord[0] = dati->coord[0];
        this_node->split.coord[1] = dati->coord[1];
        return this_node;
    }

    if (size == 2) { // onle 2 node, one null leaf and build another one node

        this_node->split.coord[0] = dati->coord[0];
        this_node->split.coord[1] = dati->coord[1];

        this_node->left = NULL;
        this_node->right = buildtree_master(++dati, 2, myaxis, --size, ++dept, max_dept);
        return this_node;
    }

    unsigned int median = ((size % 2) == 0) ? size / 2 : size / 2 + 1;

    if (size > 12) {
        TYPE median_value = find_kth(dati, median, myaxis, size);
        // median = three_way_partition(dati, median_value, size, myaxis);
    } else {
        (*sortin_network[size - 1])(dati, myaxis);
    }
    --median;

    this_node->split.coord[0] = dati[median].coord[0];
    this_node->split.coord[1] = dati[median].coord[1];

    kpoint *l_slice = dati;
    kpoint *r_slice = dati + median + 1;
    long int l_size = median;
    long int r_size = size - median - 1;

    int nthread = omp_get_num_threads();
    int myid = omp_get_thread_num();

    if (dept == max_dept) {
        if (omp_get_thread_num() == MASTER) {
            int rslave;
#pragma omp atomic capture
            rslave = ++global_slave;
            MPI_Send(&r_size, 1, MPI_LONG_INT, rslave, 0, MPI_COMM_WORLD);
            MPI_Datatype MPI_point;
            MPI_Type_contiguous(sizeof(kpoint), MPI_BYTE, &MPI_point);
            MPI_Type_commit(&MPI_point);
            MPI_Send(r_slice, r_size, MPI_point, rslave, 0, MPI_COMM_WORLD);
            knode *rnode = malloc(r_size * sizeof(knode));
#ifdef DEBUG
            sleep(1);
            for (int i = 0; i < r_size; ++i)
                printf("OMP %p [%d]r\n", &rnode[i], i);
            sleep(1);
            printf("0) -> I'm %d/%d, I will send to %d and myselft\nroffset %p\n", myid, nthread, rslave, rnode);
            fflush(stdout);
#endif
            MPI_Send(&rnode, 1, MPI_LONG_INT, rslave, 0, MPI_COMM_WORLD);
            MPI_Send(&myaxis, 1, MPI_INT, rslave, 0, MPI_COMM_WORLD);
            global_node_address = malloc(l_size * sizeof(knode));

#pragma omp parallel
            {
                int namelen;
                char processor_name[MPI_MAX_PROCESSOR_NAME];
                MPI_Get_processor_name(processor_name, &namelen);
                int np = omp_get_num_threads();
                int iam = omp_get_thread_num();
                int cpu_num = sched_getcpu();
                fprintf(stderr, "core %d -> %s | MASTER tree | th %d/%d \n", cpu_num, processor_name, iam, np);

#pragma omp single
                {
                    this_node->left = (buildtree(l_slice, 2, myaxis, l_size, NULL, NULL));
                }
            }
            MPI_Datatype MPI_node;
            MPI_Type_contiguous(sizeof(knode), MPI_BYTE, &MPI_node);
            MPI_Type_commit(&MPI_node);

            MPI_Status stat;
            MPI_Recv(rnode, r_size, MPI_node, rslave, 0, MPI_COMM_WORLD, &stat);
            this_node->right = rnode;

        } else {
            int rslave, lslave;
#pragma omp atomic capture
            rslave = ++global_slave;
#pragma omp atomic capture
            lslave = ++global_slave;

            MPI_Send(&r_size, 1, MPI_LONG_INT, rslave, 0, MPI_COMM_WORLD);
            MPI_Send(&l_size, 1, MPI_LONG_INT, lslave, 0, MPI_COMM_WORLD);
            MPI_Datatype MPI_point;
            MPI_Type_contiguous(sizeof(kpoint), MPI_BYTE, &MPI_point);
            MPI_Type_commit(&MPI_point);
            // create type to send point
            MPI_Send(r_slice, r_size, MPI_point, rslave, 0, MPI_COMM_WORLD);
            MPI_Send(l_slice, l_size, MPI_point, lslave, 0, MPI_COMM_WORLD);

            knode *lnode = malloc(l_size * sizeof(knode));
            knode *rnode = malloc(r_size * sizeof(knode));
#ifdef DEBUG
            for (int i = 0; i < r_size; ++i)
                printf("OMP %p [%d]r\n", &rnode[i], i);
            for (int i = 0; i < l_size; ++i)
                printf("OMP %p [%d]l\n", &lnode[i], i);
#endif
            MPI_Send(&lnode, 1, MPI_LONG_INT, lslave, 0, MPI_COMM_WORLD);
            MPI_Send(&rnode, 1, MPI_LONG_INT, rslave, 0, MPI_COMM_WORLD);
            MPI_Send(&myaxis, 1, MPI_INT, rslave, 0, MPI_COMM_WORLD);
            MPI_Send(&myaxis, 1, MPI_INT, lslave, 0, MPI_COMM_WORLD);
            MPI_Status stat;
            // create type to receive node
            MPI_Datatype MPI_node;
            MPI_Type_contiguous(sizeof(knode), MPI_BYTE, &MPI_node);
            MPI_Type_commit(&MPI_node);
            MPI_Recv(rnode, r_size, MPI_node, rslave, 0, MPI_COMM_WORLD, &stat);
            MPI_Recv(lnode, l_size, MPI_node, lslave, 0, MPI_COMM_WORLD, &stat);
            this_node->right = rnode;
            this_node->left = lnode;
#ifdef DEBUG
            printf("0) I'm %d/%d, I will send to %d and %d\n roffset %p\n loffset %p\n", myid, nthread, rslave, lslave, rnode, lnode);
            fflush(stdout);
#endif
        }

    } else {
#pragma omp task shared(dati) firstprivate(myaxis, median, size, dept, max_dept)
        this_node->left = buildtree_master(l_slice, 2, myaxis, l_size, ++dept, max_dept);
#pragma omp task shared(dati) firstprivate(myaxis, median, size, dept, max_dept)
        this_node->right = buildtree_master(r_slice, 2, myaxis, r_size, ++dept, max_dept);
    }
    return this_node;
}

int main(int argc, char *argv[]) {

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided < MPI_THREAD_MULTIPLE) {
        printf("Sta cosa del thread multiple non si puo fa ! \n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (argc < 4) {
        if (rank == MASTER)
            printf("Not enought args\n!");
        MPI_Finalize();
        return -1;
    }
    int namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name, &namelen);
    omp_set_num_threads(atoi(argv[3]));
    if (rank == MASTER) {

        unsigned int p_number = atol(argv[1]);
        kpoint *my_data = malloc(sizeof(kpoint) * p_number);
        init_point(my_data, p_number);
        knode *mynode = NULL;
        omp_set_nested(1);
        omp_set_num_threads(size / 2);
#pragma omp parallel
        {
            omp_set_num_threads(atoi(argv[3]));
            int np = omp_get_num_threads();
            int iam = omp_get_thread_num();
            int cpu_num = sched_getcpu();
            fprintf(stderr, "core %d -> %s |MASTER rank %d/%d | th %d/%d |\n", cpu_num, processor_name, rank, size, iam, np);

#pragma omp single
            {
                mynode = buildtree_master(my_data, 2, -1, p_number, 0, log2l(size) - 1);
            }
        }

        if (atoi(argv[2]) == 1) {
            printf("\ngraph test {\n");
            view_tree(mynode);
            printf("}\n");
        }
        MPI_Finalize();
        return 0;
    } else {

        MPI_Status stat;
        long int data_size;

        MPI_Recv(&data_size, 1, MPI_LONG_INT, MASTER, 0, MPI_COMM_WORLD, &stat);
        kpoint *data = malloc(data_size * sizeof(kpoint));
        MPI_Datatype MPI_point;
        MPI_Type_contiguous(sizeof(kpoint), MPI_BYTE, &MPI_point);
        MPI_Type_commit(&MPI_point);
        MPI_Recv(data, data_size, MPI_point, MASTER, 0, MPI_COMM_WORLD, &stat);
        knode *memory_offset = NULL;
        MPI_Recv(&memory_offset, 1, MPI_LONG_INT, MASTER, 0, MPI_COMM_WORLD, &stat);
        int axis;
        MPI_Recv(&axis, 1, MPI_INT, MASTER, 0, MPI_COMM_WORLD, &stat);
        /* #ifdef DEBUG */
        /*     printf("1) I'm %d and my size is %d\n myoffset %p\n",rank,data_size,memory_offset); */
        /*     fflush(stdout); */
        /* #endif */
        global_node_address = malloc(data_size * sizeof(knode));
        knode *mynode = NULL;

        omp_set_num_threads(atoi(argv[3]));
#pragma omp parallel
        {
            int np = omp_get_num_threads();
            int iam = omp_get_thread_num();
            int cpu_num = sched_getcpu();

            fprintf(stderr, "core %d -> %s |rank %d/%d | th %d/%d |\n", cpu_num, processor_name, rank, size, iam, np);

#pragma omp single
            {
                mynode = buildtree(data, NDIM, axis, data_size, global_node_address, memory_offset);
            }
        }
        MPI_Datatype MPI_node;
        MPI_Type_contiguous(sizeof(knode), MPI_BYTE, &MPI_node);
        MPI_Type_commit(&MPI_node);

        MPI_Send(mynode, data_size, MPI_node, MASTER, 0, MPI_COMM_WORLD);

        MPI_Finalize();
    }
}
