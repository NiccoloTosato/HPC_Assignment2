#define _GNU_SOURCE //needed to get binding information
#include <sched.h>

#include <mpi.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#define NDIM 2
#define TYPE double
#define MIN_SIZE 10
typedef struct knode knode;
typedef struct kpoint kpoint;



struct kpoint {
  TYPE coord[NDIM];
};

struct knode {
  kpoint split;
  knode* left;
  knode *right;
  int axis;

};



knode* global_address=NULL;
knode* base_address=NULL;
knode* memory_offset=NULL;
int slave_number=1;

void  init_point(kpoint* dati,int n){  
  //short unsigned int state[3];
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
  
  
  TYPE pivot = array[high].coord[axis];
  //low element
  int i = (low - 1);

  //partition
  for (int j = low; j < high; j++) {
    if (array[j].coord[axis] <= pivot) {      
      i++;
      swap(&array[i], &array[j]);
    }
  }

  swap(&array[i + 1], &array[high]);  
  // return the partition point
  return (i + 1);
}

void quickSort(kpoint* array, int low, int high,int axis) {

  if (low < high) {
    
    int pi = partition(array, low, high, axis);
    quickSort(array, low, pi - 1, axis);
    quickSort(array, pi + 1, high, axis);
  }
}


void view_tree(knode* dati){
  if ( (dati->right==NULL) && (dati->left==NULL) ){
    return;
  }
  
  if (dati->right==NULL) {
    kpoint point=(dati->split);
    knode *left= dati->left;

    kpoint l_point = (left->split);
    printf(" \"%f,%f\" -- \"%f,%f\" [label=%c] \n",point.coord[0],point.coord[1],l_point.coord[0],l_point.coord[1],dati->axis+88);

    view_tree(left);
    return;
  }
  
  if(dati->left==NULL) {
    kpoint point=(dati->split);
    knode *right= dati->right;
    kpoint r_point = (right->split);
    printf(" \"%f,%f\" -- \"%f,%f\" [label=%c]\n",point.coord[0],point.coord[1],r_point.coord[0],r_point.coord[1],dati->axis+88);

    view_tree(right);
    return;
  }
  
  kpoint point=(dati->split);
  knode *right= dati->right;
  knode *left= dati->left;
  kpoint r_point = (right->split);
  kpoint l_point = (left->split);

  printf(" \"%f,%f\" -- \"%f,%f\" [label=%c] \n",point.coord[0],point.coord[1],r_point.coord[0],r_point.coord[1],dati->axis+88);
  printf(" \"%f,%f\" -- \"%f,%f\" [label=%c] \n",point.coord[0],point.coord[1],l_point.coord[0],l_point.coord[1],dati->axis+88);
  
  view_tree(left);
  view_tree(right);
  return;
  

}

void  view_node(knode* dati) {
  kpoint point=(dati->split);
  printf("(%f,%f), axis %d\n",point.coord[0],point.coord[1],1);
}

int  three_way_partition(kpoint* dati,TYPE midvalue, int npoint,int axis){
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

  return j;
}


TYPE find_kth(kpoint* dati,int k,int axis,int npoints) {


  TYPE pivot=dati[0].coord[axis];

  int split_point=three_way_partition(dati, pivot, npoints, axis);

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



  knode* this_node=NULL;
  #pragma omp atomic capture
    this_node=++global_address;
  //set axis
  int myaxis = (axis+1) % ndim;    
  this_node->axis=myaxis;



  if( size > 2) {
    int median=size/2;//=((size %2) ==0)? size/2:size/2+1;
    
    if( size > MIN_SIZE ) {
      TYPE median_value=find_kth(dati, median, myaxis, size);
      median=three_way_partition(dati, median_value, size, myaxis);
      --median;
    } else {
      quickSort(dati, 0, size-1, myaxis);
    }

    ////////////////////STANDARD WAY/////////////
    this_node->split.coord[0]=dati[median].coord[0];
    this_node->split.coord[1]=dati[median].coord[1];
        
    //////////////////ALTERNATIVE WAY///////////////////
    /* don't move data,save address */
    /* this_node->split=dati+median ; */
    ////////////////////////////////////////////////////


    #pragma omp task shared(dati) firstprivate(myaxis, median,size)
    {
      this_node->left =  (buildtree(dati, 2, myaxis, median)-base_address)*sizeof(knode)+(long int) memory_offset;
      printf("MPI left %p\n",this_node->left);
    }
    #pragma omp task shared(dati) firstprivate(myaxis, median,size)
    {
      this_node->right = (buildtree(dati+median+1, 2, myaxis,size-median-1)-base_address)*sizeof(knode) + (long int)memory_offset;
      printf("MPI right %p\n",this_node->right);
    }
    return this_node;
  } else if (size == 1) { 

    
    //////////////////ALTERNATIVE WAY///////////////////
    /* don't move data,save only address */
    /* this_node->split=*dati; */
    ////////////////////////////////////////////////////
    /////////////////STANDARD WAY///////////////////////
    this_node->split.coord[0]=dati->coord[0];
    this_node->split.coord[1]=dati->coord[1];
    ///////////////////////////////////////////////////

    this_node->left = NULL;
    this_node->right = NULL;
    
    return this_node;
    
  } else   if (size ==2) {

    /////////////////STANDARD WAY///////////////////////
    this_node->split.coord[0]=dati->coord[0];
    this_node->split.coord[1]=dati->coord[1];
    ////////////////////////////////////////////////////
    //////////////////ALTERNATIVE WAY///////////////////
    /* don't move data,save only address */
    /* this_node->split=*dati; */
    ////////////////////////////////////////////////////
    
    this_node->left = NULL;
    this_node->right = (buildtree(dati+1, 2, myaxis, 1)-base_address)*sizeof(knode) + (long int) memory_offset;
    printf("MPI left %p\n",this_node->right);
    return this_node;
    
  }

}


knode*   buildtree_master(kpoint* dati,int ndim,int axis,int size,int dept,int worker){
  
  knode* this_node=malloc(sizeof(knode));
  kpoint* my_split=(kpoint*) malloc(sizeof(kpoint));
  int myaxis = (axis+1) % ndim;
  this_node->axis=myaxis;
      
  int median=size/2;
  
  if(size>MIN_SIZE){
    TYPE median_value=find_kth(dati, median, myaxis, size);
    median=three_way_partition(dati, median_value, size, myaxis);
    --median;
  } else {
    quickSort(dati, 0, size-1, myaxis);
  }

  //save data in split
  my_split->coord[0]=dati[median].coord[0];
  my_split->coord[1]=dati[median].coord[1];
  this_node->split=*my_split;

  /* instead don't move data,save address */
  /* this_node->split=dati+median ; */
  
  if(pow(2,dept+1)==worker-1) {
    int my_slave_r;
    int my_slave_l;
    
    #pragma omp atomic capture  
    my_slave_r=slave_number++;
    #pragma omp atomic capture  
    my_slave_l=slave_number++;
    int r_size,l_size;
    kpoint* r_slice;
    kpoint* l_slice;
    int my_tag=omp_get_thread_num();
    r_size=size-median-1;
    l_size=median;
    
    r_slice=dati+median+1;
    l_slice=dati;
    
    
    MPI_Send(&r_size, 1, MPI_INT, my_slave_r, 0, MPI_COMM_WORLD);
    MPI_Send(&l_size, 1, MPI_INT, my_slave_l, 0, MPI_COMM_WORLD);

    MPI_Send(r_slice, r_size*sizeof(kpoint), MPI_BYTE, my_slave_r, 0, MPI_COMM_WORLD);
    MPI_Send(l_slice, l_size*sizeof(kpoint), MPI_BYTE, my_slave_l, 0, MPI_COMM_WORLD);
    
    knode* rnode=calloc(sizeof(knode),r_size);
    knode* lnode=calloc(sizeof(knode),l_size);
    printf("rlSIZE %d %d\n",r_size,l_size);
    printf("this_node  at %p, th %d\n",this_node,omp_get_thread_num());
    for(int i=0;i<r_size;++i)
      printf("rnode at %p, th %d\n",&rnode[i],omp_get_thread_num());
    for(int i=0;i<l_size;++i)
      printf("lnode at %p, th %d\n",&lnode[i],omp_get_thread_num());

    ////////////////////////////SEND ADDRESS PART//////////////////////////
    //long int rnode_int=(long int) rnode;
    //long int lnode_int=(long int) lnode;
    MPI_Send(&rnode, 1, MPI_LONG_INT, my_slave_r, 0, MPI_COMM_WORLD);
    MPI_Send(&lnode, 1, MPI_LONG_INT, my_slave_l, 0, MPI_COMM_WORLD);
    ///////////////////////////////////////////////////////////////////////
    
    //printf("I'm %d, right %p | left %p\n",omp_get_thread_num(),rnode,lnode);
    MPI_Request req[2];
    MPI_Status stat[2];


      MPI_Irecv(rnode, sizeof(knode)*r_size, MPI_BYTE, my_slave_r, MPI_ANY_TAG, MPI_COMM_WORLD, &req[0]);
      MPI_Irecv(lnode, sizeof(knode)*l_size, MPI_BYTE, my_slave_l, MPI_ANY_TAG, MPI_COMM_WORLD, &req[1]);

      MPI_Wait(&req[0], &stat[0]);
      MPI_Wait(&req[1], &stat[1]);
      //printf("I'm %d, received all\n",omp_get_thread_num());

    //printf("I'reaced %d dept, STOPHERE ,I'm %d\nI send to %d and %d\n",dept,omp_get_thread_num(),my_slave_r,my_slave_l);
    this_node->left=lnode;
    this_node->right=rnode;


  } else {
#pragma omp task shared(dati) firstprivate(myaxis, median,size,worker,dept)
{
  //printf("I'm %d | continue build dept %d\n",omp_get_thread_num(),dept);
  this_node->left = buildtree_master(dati, 2, myaxis, median,dept+1,worker);
 }
#pragma omp task shared(dati) firstprivate(myaxis, median,size,worker,dept)
 {
   //printf("I'm %d | continue build dept %d\n",omp_get_thread_num(),dept);
   this_node->right = buildtree_master(dati+median+1, 2, myaxis,size-median-1,dept+1,worker);
 }
  
  }
  return this_node;
}



#define MASTER 0
int   main(int argc, char* argv[]) {

  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

  if(provided < MPI_THREAD_MULTIPLE)
    {
      printf("The threading support level is lesser than that demanded.\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  int size,rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  


  if(rank==MASTER) {
    if (argc<1){
      printf("Not enought args\n!");
      return -1;
    }
    int p_number=atoi(argv[1]);
    kpoint* my_data=malloc(sizeof(kpoint)*p_number);
    omp_set_num_threads(8);
          init_point(my_data,p_number);
/* #pragma omp parallel */
/*     { */

/*     } */
    //printf("fine init\n");
  
    //knode* mynode=malloc(sizeof(knode)*p_number);
    knode* mynode=malloc(sizeof(knode));
    global_address=mynode-1;
    omp_set_num_threads(size-1);


    int namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name, &namelen);

    omp_set_num_threads(size);
    #pragma omp parallel
    {
      int np,iam;

      np = omp_get_num_threads();
      iam = omp_get_thread_num();
      int cpu_num = sched_getcpu();
      printf("rank %d/%d | th %d/%d | core %d |node %s\n",rank,size,iam,np,cpu_num,processor_name);

      #pragma omp single
      {
	mynode=buildtree_master(my_data, 2, -1, p_number,0,size);
      }
    }


    if(atoi(argv[2])==1){
      printf("\ngraph test {\n");
      view_tree(mynode);
      printf("}\n");
    }

    
  } else {
    //SLAVE PART !
    MPI_Status stat;
    int data_size;
    MPI_Recv(&data_size, 1, MPI_INT,MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);

    int namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(processor_name, &namelen);
    kpoint* my_data=calloc(data_size,sizeof(kpoint));
    knode* mynode=calloc(data_size,sizeof(knode));
    global_address=mynode-1;
    base_address=mynode;

    //int my_master_tag=stat.MPI_TAG;

    MPI_Recv(my_data, data_size*sizeof(kpoint),MPI_BYTE, MASTER, MPI_ANY_TAG,MPI_COMM_WORLD, &stat);
    //printf("RANK %d | data_size %d\n",rank,data_size);

    //long int memory_offset_int;
    MPI_Recv(&memory_offset, 1, MPI_LONG_INT, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
    //memory_offset=(knode*)memory_offset_int;
    //printf("my rank %d,my offset %p, BASEADDRESS %p\n",rank,memory_offset,mynode);
    int num_th=atoi(argv[3]);
    omp_set_num_threads(num_th);
    #pragma omp parallel
    {
      int np,iam;

      np = omp_get_num_threads();
      iam = omp_get_thread_num();
      int cpu_num = sched_getcpu();
      printf("rank %d/%d | th %d/%d | core %d |node %s\n",rank,size,iam,np,cpu_num,processor_name);

      #pragma omp single
      {
    
	mynode=buildtree(my_data, 2, -1, data_size);
      }
    }
    MPI_Send(mynode, data_size*sizeof(knode), MPI_BYTE, MASTER, 0, MPI_COMM_WORLD);

  }
  MPI_Finalize();
}


