#include "/usr/lib/x86_64-linux-gnu/openmpi/include/mpi.h"
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <unistd.h>

#define NDIM 2
#define NPOINT 12
typedef struct knode knode;
typedef struct kpoint kpoint;
#define TYPE double
struct kpoint {
  TYPE coord[NDIM];
};

struct knode {
  kpoint split;
  knode* left;
  knode *right;

  int axis;

};



void init_point(kpoint* dati,int n){
  //#pragma omp parallel for
  for(int i=0;i<n;++i){
    dati[i].coord[0]=drand48();
    dati[i].coord[1]=drand48();      
  }
}


void swap(kpoint* restrict a, kpoint* restrict b) {
  kpoint t = *a;
  a->coord[0]=b->coord[0];
  a->coord[1]=b->coord[1];
  
  b->coord[0]=t.coord[0];
  b->coord[1]=t.coord[1];
}

int partition(kpoint* array, int low, int high,int axis) {
  
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

void quickSort(kpoint* array, int low, int high,int axis) {

  if (low < high) {
    
    // find the pivot element such that
    // elements smaller than pivot are on left of pivot
    // elements greater than pivot are on right of pivot
    int pi = partition(array, low, high, axis);
    //printf("Hi, i'm %d/%d in sorting\n",omp_get_thread_num(),omp_get_num_threads());    
    // recursive call on the left of pivot

    quickSort(array, low, pi - 1, axis);
    
    // recursive call on the right of pivot

    quickSort(array, pi + 1, high, axis);
  }
}


void view_tree(knode* dati){
  if ( (dati->right==NULL) && (dati->left==NULL) ){
    return;
  }
  
  if (dati->right==NULL) {
    kpoint point=dati->split;
    knode *left= dati->left;
    kpoint l_point = left->split;
    printf(" \"%f,%f\" -- \"%f,%f\" [label=%c] \n",point.coord[0],point.coord[1],l_point.coord[0],l_point.coord[1],dati->axis+88);

    view_tree(left);
    return;
  }
  
  if(dati->left==NULL) {
    kpoint point=dati->split;
    knode *right= dati->right;
    kpoint r_point = right->split;
    printf(" \"%f,%f\" -- \"%f,%f\" [label=%c]\n",point.coord[0],point.coord[1],r_point.coord[0],r_point.coord[1],dati->axis+88);

    view_tree(right);
    return;
  }
  
  kpoint point=dati->split;
  knode *right= dati->right;
  knode *left= dati->left;
  
  kpoint r_point = right->split;
  kpoint l_point = left->split;

  printf(" \"%f,%f\" -- \"%f,%f\" [label=%c] \n",point.coord[0],point.coord[1],r_point.coord[0],r_point.coord[1],dati->axis+88);
  printf(" \"%f,%f\" -- \"%f,%f\" [label=%c] \n",point.coord[0],point.coord[1],l_point.coord[0],l_point.coord[1],dati->axis+88);
  
  view_tree(left);
  view_tree(right);
  return;
  

}

void view_node(knode* dati) {
  kpoint point=dati->split;
  printf("(%f,%f), axis %d\n",point.coord[0],point.coord[1],1);
}

int three_way_partition(kpoint* dati,TYPE midvalue, int npoint,int axis){
  int i=0;
  int j=0;
  int k=npoint -1;
  while(j<=k){
    if(dati[j].coord[axis] < midvalue) {
      swap(dati+i,dati+j);
      ++i;
      ++j;
    } else if(dati[j].coord[axis] > midvalue){
      swap(dati+j,dati+k);
      --k;
    } else {
      ++j;
    }

  }
  //printf("j %d, k %d\n",j,k);
  return j;
}

TYPE find_mean(kpoint* dati,int npoint,int axis) {
  TYPE min=dati[0].coord[axis];
  TYPE max=dati[0].coord[axis];
  for(int i=0;i<npoint;++i){
    if(dati[i].coord[axis]> max){
      max=dati[i].coord[axis];
      continue;
    }
    if(dati[i].coord[axis]< min){
      min=dati[i].coord[axis];
    }
  }
  

  return (min+max)/2;
}

int find_pivot(kpoint* dati,TYPE mean,TYPE npoint,int axis){
  for(int i=0;i<npoint;++i){
    if(dati[i].coord[axis] == mean)
      return (i);
  }
  return -1;
}

TYPE find_kth(kpoint* dati,int k,int axis,int npoints){

  //printf("\nk %d\n",k);
  TYPE pivot=dati[0].coord[axis];

  int split_point=three_way_partition(dati, pivot, npoints, axis);
  //printf("sp %d ",split_point);
  //split_point=find_pivot(dati, pivot, npoints, axis);
  //printf("sp2 %d \n",split_point);
  kpoint* l=dati;
  kpoint* r=dati+split_point;

  int size_r=split_point-1;
  int size_l=npoints-split_point;
  
  if(k==size_l+1) return pivot;
  if(k<size_l+1) return find_kth(r, k, axis, size_l);
  if(k>size_l+1) return find_kth(l, k-(size_l+1), axis, size_r);
  return -1;
}

knode* buildtree(kpoint* dati,int ndim,int axis,int size){
  //printf("INIT SIZE %d\n",size);
  int myaxis = (axis+1) % ndim;
  knode* this_node= (knode*) malloc(sizeof(knode));
  this_node->axis=myaxis;
  kpoint* my_split=(kpoint*) malloc(sizeof(kpoint));

  if (size == 1) {
    this_node->left = NULL;
    this_node->right = NULL;
    my_split->coord[0]=dati->coord[0];
    my_split->coord[1]=dati->coord[1];
    this_node->split=*my_split;

    return this_node;
  }
  
  if (size ==2) {

    my_split->coord[0]=dati->coord[0];
    my_split->coord[1]=dati->coord[1];
    this_node->split=*my_split;
    this_node->left = NULL;

    this_node->right = buildtree(dati+1, 2, myaxis, 1);
    
    return this_node;
    
  }
  /* double mean=find_extreme(dati, size, myaxis); */
  /* three_way_partition(dati, mean, size, myaxis); */
  int median=size/2;//((size %2) ==0) ? size/2 : size/2+1;
 printf("%d median\n",median);
  if(size>12){


    TYPE median_value=find_kth(dati, median, myaxis, size);
    median=three_way_partition(dati, median_value, size, myaxis);
    //printf("\nmedian %d, median value%f,size %d\n",median,median_value,size);

    //view_point(dati, size);
    //printf("\n\n");
    --median;
  } else {
    //printf("quicksort call\n");
    quickSort(dati, 0, size-1, myaxis);
  }

  //printf("MEDIAN %d,SIZE %d\n",median,size);

  if(size==0){
    
    sleep(99999);
    printf("STOP\n");
  }
  my_split->coord[0]=dati[median].coord[0];
  my_split->coord[1]=dati[median].coord[1];

  this_node->split=*my_split;

  
  //printf("I'm %d/%d",omp_get_thread_num(),omp_get_num_threads());
  // 11 elementi, mediano 6, size/2 = 5
  #pragma omp task shared(dati) firstprivate(myaxis, median,size)
  this_node->left = buildtree(dati, 2, myaxis, median);
  #pragma omp task shared(dati) firstprivate(myaxis, median,size)
  this_node->right = buildtree(dati+median+1, 2, myaxis,size-median-1);
  return this_node;
}



int main(int argc, char* argv[]) {
  if (argc<1){
    printf("Not enought args\n!");
    return -1;
  }
  int p_number=atoi(argv[1]);
  kpoint* my_data=malloc(sizeof(kpoint)*p_number);
  init_point(my_data,p_number);
  //printf("fine init\n");

  knode *mynode;
  #pragma omp parallel
  {
    #pragma omp single
      {
	mynode=buildtree(my_data, 2, -1, p_number);
      }
  }
  /* mynode=buildtree(my_data, 2, -1, p_number, 0); */
  if(atoi(argv[2])==1){
  printf("\ngraph test {\n");
  view_tree(mynode);
  printf("}\n");
    }
}

