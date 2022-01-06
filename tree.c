#include "/usr/lib/x86_64-linux-gnu/openmpi/include/mpi.h"
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>


#define NDIM 2
typedef struct knode knode;
typedef struct kpoint kpoint;

struct kpoint {
  int coord[NDIM];
};

struct knode {
  int axis;
  kpoint split;
  knode* left;
  knode *right;
};


void view_point(kpoint* dati,int n){
  for(int i=0;i<n;i++){
    printf("(%d,%d), ",dati[i].coord[0],dati[i].coord[1]);
  }
  printf("\n");
}

// function to swap elements
void swap(kpoint* a, kpoint* b) {
  kpoint t = *a;
  a->coord[0]=b->coord[0];
  a->coord[1]=b->coord[1];
  
  b->coord[0]=t.coord[0];
  b->coord[1]=t.coord[1];
}

// function to find the partition position
int partition(kpoint* array, int low, int high,int axis) {
  
  // select the rightmost element as pivot
  int pivot = array[high].coord[axis];
  
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
#pragma omp task shared(array) firstprivate(low,pi,axis)
    quickSort(array, low, pi - 1, axis);
    
    // recursive call on the right of pivot
#pragma omp task shared(array) firstprivate(high,pi,axis)
    quickSort(array, pi + 1, high, axis);
  }
}


knode* buildtree(kpoint* dati,int ndim,int axis,int size,int dept){
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
    //printf("\nsize==1, ");
    //printf("Buildtree before split:\n");
    //view_point(dati, size);

    return this_node;
  }
  if (size ==2) {
    
    my_split->coord[0]=dati->coord[0];
    my_split->coord[1]=dati->coord[1];
    this_node->split=*my_split;
    this_node->left = NULL;
    this_node->right = buildtree(dati+1, 2, myaxis, 1,dept+1);
    //printf("\nsize==2, ");
    //printf("Buildtree before split:\n");
    //view_point(dati, size);
    
    return this_node;
    
  }
  if(dept==1){
    omp_set_num_threads(8);

  }
  if(dept==2)
    omp_set_num_threads(4);
  if(dept>2)
    omp_set_num_threads(1);
#pragma omp parallel
  {
    #pragma omp single
    {    
      //printf("8 th,%d\n",omp_get_num_threads());
      quickSort(dati,0,size-1,myaxis); // order all eleme
    }
}
  //printf("\nBuildtree before split:\n");
  //view_point(dati, size);
  int median= size/2; //(size % 2 == 0 ) ? size / 2 -1: size / 2 ;
  //printf("Splitpoint: (%d,%d),median %d\n",dati[median].coord[0],dati[median].coord[1],median);


  my_split->coord[0]=dati[median].coord[0];
  my_split->coord[1]=dati[median].coord[1];

  this_node->split=*my_split;

  
  //printf("I'm %d/%d",omp_get_thread_num(),omp_get_num_threads());
  // 11 elementi, mediano 6, size/2 = 5
#pragma omp task shared(dati) firstprivate(myaxis, median,size,dept)
  this_node->left = buildtree(dati, 2, myaxis, median ,dept+1);
#pragma omp task shared(dati) firstprivate(myaxis, median,size,dept)
  this_node->right = buildtree(dati+median+1, 2, myaxis,size-median-1 ,dept+1);
  return this_node;
}


void init_point(kpoint* dati,int n){
  for(int i=0;i<n;++i){
    dati[i].coord[0]=rand() % 100;
    dati[i].coord[1]=rand() % 100;      
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
    printf(" \"%d,%d\" -- \"%d,%d\" [label=%d] \n",point.coord[0],point.coord[1],l_point.coord[0],l_point.coord[1],dati->axis);

    view_tree(left);
    return;
  }
  
  if(dati->left==NULL) {
    kpoint point=dati->split;
    knode *right= dati->right;
    kpoint r_point = right->split;
    printf(" \"%d,%d\" -- \"%d,%d\" [label=%d]\n",point.coord[0],point.coord[1],r_point.coord[0],r_point.coord[1],dati->axis);

    view_tree(right);
    return;
  }
  
  kpoint point=dati->split;
  knode *right= dati->right;
  knode *left= dati->left;
  
  kpoint r_point = right->split;
  kpoint l_point = left->split;

  printf(" \"%d,%d\" -- \"%d,%d\" [label=%c] \n",point.coord[0],point.coord[1],r_point.coord[0],r_point.coord[1],dati->axis+88);
  printf(" \"%d,%d\" -- \"%d,%d\" [label=%c] \n",point.coord[0],point.coord[1],l_point.coord[0],l_point.coord[1],dati->axis+88);
  
  view_tree(left);
  view_tree(right);
  return;
  

}

void view_node(knode* dati) {
  kpoint point=dati->split;
  printf("(%d,%d), axis %d\n",point.coord[0],point.coord[1],1);
}
int main(int argc, char* argv[]) {
  srand(10);
  int p_number=atoi(argv[1]);
  kpoint* my_data=malloc(sizeof(kpoint)*p_number);
  init_point(my_data,p_number);
  //quickSort(my_data, 0, p_number-1, 1);
  //view_point(my_data,NPOINT);
  knode* mynode;


  #pragma omp parallel
  {
    omp_set_nested(1);
    //printf("omp nested: %d\n",omp_get_nested());
    #pragma omp single
      {
	mynode = buildtree(my_data, 2, -1, p_number,1);
      }
  }
  //view_node(mynode);
  /* printf("graph test {\n"); */
  /* view_tree(mynode); */
  /* printf("}\n"); */
}
// 1.46 1.45 1.47 | .38 .38
