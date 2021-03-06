
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define NDIM 2
typedef struct knode knode;
typedef struct kpoint kpoint;
#define TYPE double

struct kpoint {
    TYPE coord[NDIM];
};

struct knode {
    kpoint split;
    knode *left;
    knode *right;
    unsigned short int axis;
};
knode *global_node_address = NULL;

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

int partition(kpoint *array, int low, int high, short int axis) {

    // select the rightmost element as pivot
    TYPE pivot = array[high].coord[axis];

    // pointer for greater element
    int i = (low - 1);

    // traverse each element of the array
    // compare them with the pivot
    for (int j = low; j < high; j++) {
        if (array[j].coord[axis] <= pivot) {

            // if element smaller than pivot is found
            // swap it with the greater element pointed by i
            i++;

            // swap element at i with element at j
            swap(&array[i], &array[j]);
        }
    }

    // swap the pivot element with the greater element at i
    swap(&array[i + 1], &array[high]);

    // return the partition point
    return (i + 1);
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

TYPE find_kth(kpoint *dati, unsigned int k, unsigned short int axis, unsigned int npoints) {

    // printf("\nk %d\n",k);
    TYPE pivot = dati[0].coord[axis];

    unsigned int split_point = three_way_partition(dati, pivot, npoints, axis);
    // printf("sp %d ",split_point);
    // split_point=find_pivot(dati, pivot, npoints, axis);
    // printf("sp2 %d \n",split_point);
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

knode *buildtree(kpoint *dati, int ndim, short int axis, unsigned int size) {
    int myaxis = (axis + 1) % ndim;
    knode *this_node = NULL;
#pragma omp atomic capture
    this_node = global_node_address++;
    this_node->axis = myaxis;
    if (size == 1) {
        this_node->left = NULL;
        this_node->right = NULL;
        this_node->split.coord[0] = dati->coord[0];
        this_node->split.coord[1] = dati->coord[1];
        return this_node;
    }

    if (size == 2) {

        this_node->split.coord[0] = dati->coord[0];
        this_node->split.coord[1] = dati->coord[1];
        this_node->right = buildtree(dati + 1, 2, myaxis, 1);
        this_node->left = NULL;

        return this_node;
    }

    unsigned int median = ((size % 2) == 0) ? size / 2 : size / 2 + 1;

    TYPE median_value = find_kth(dati, median, myaxis, size);
   

    --median;

    this_node->split.coord[0] = dati[median].coord[0];
    this_node->split.coord[1] = dati[median].coord[1];

    kpoint *l_slice = dati;
    kpoint *r_slice = dati + median + 1;
    unsigned int l_size = median;
    unsigned int r_size = size - median - 1;

#pragma omp task shared(dati) // firstprivate(myaxis, median, size)
    this_node->left = (buildtree(l_slice, 2, myaxis, l_size));
#pragma omp task shared(dati) // firstprivate(myaxis, median, size)
    this_node->right = (buildtree(r_slice, 2, myaxis, r_size));
    return this_node;
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("Not enought args\n!");
        return -1;
    }
    omp_set_num_threads(atoi(argv[3]));

    unsigned int p_number = atoi(argv[1]);

    kpoint *my_data = malloc(sizeof(kpoint) * p_number);
    init_point(my_data, p_number);
    global_node_address = malloc(p_number * sizeof(knode));

    knode *mynode;
#pragma omp parallel
    {

#pragma omp single
        {
            mynode = buildtree(my_data, 2, -1, p_number);
        }
    }

    if (atoi(argv[2]) == 1) {
        printf("\ngraph test {\n");
        view_tree(mynode);
        printf("}\n");
    }
}
