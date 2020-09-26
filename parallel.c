/*
* FILE: parallel.c
* THMMY, 7th semester, Parallel and Distributed Systems: 4th assignment
* Parallel Implementation of the Reverse Cuthill Mckee Algorithm
* Author:
*   Moustaklis Apostolos, 9127, amoustakl@auth.gr
* It will create a N X N Sparse Matrix
* And transform it into a band matrix form with a small bandwidth.
*/

#include <stdio.h>
#include <stdlib.h>
#include "queue.h"
#include <time.h>
#include <sys/time.h>
#include <omp.h>



#define N 10000
#define SPARSITY_LIMIT 0.8
#define QUEUE_SIZE 10000
#define MAX_THREADS 8
#define SAVE 0


//Function Decleration
void MatrixInitialization(int * matrix , int size , double sparsityLimit);
void degreeCalculator(int *matrix , int * degreeArray);
void ReverseCuthilMckee(int *matrix , int *degreeArray , int * R );
void swapElement(int *one, int  *two);
void sortByDegree(int *neighbors ,int  *degreeArray , int neighborsCounter);
int* matrixReorder(int* matrix, int* R, int size);
int findMinIdxParallel(int* degrees, int* notVisited, int size, int threads_num);
void mergeSort_degrees_indexes(int * arr,int * idx, int l, int r);
void merge(int * arr,int * idx, int l, int m, int r);
void outputWrite(int* matrix, int rows, int col, char* file_path);

//Calculating the chunksize for the threads
int CHUNKSIZE = N / MAX_THREADS;

int main(int  argc, const char* argv[]){

  int *  matrix = (int * )calloc(N*N,sizeof(int));
  int * degreeArray = (int *)calloc(N,sizeof(int));
  int  * R = (int *)calloc(N,sizeof(int));

  if(matrix == NULL || degreeArray == NULL || R == NULL){
    printf("Error at memory Initialization \n");
    exit(0);
  }

MatrixInitialization(matrix, N , SPARSITY_LIMIT);
struct timeval start, end;


gettimeofday(&start, NULL);
degreeCalculator(matrix , degreeArray) ;
ReverseCuthilMckee(matrix , degreeArray , R);
gettimeofday(&end, NULL);

double time = ((double)((end.tv_sec*1e6 + end.tv_usec) - (start.tv_sec*1e6 + start.tv_usec)))*1e-6;
printf("Execution Time: %lf sec\n", time);

int save = SAVE;
if ( save){
   int* reorderedMatrix = matrixReorder( matrix, R, N) ;
   outputWrite(reorderedMatrix, N, N, "outputMatrix.txt");
}


};


//Matrix Initialization function
void MatrixInitialization (int * matrix , int size , double sparsityLimit  ){

//sparsityLimit the percentage of the zeros
  int nonZeros = (size*size) - (size*size*sparsityLimit);
  int sum = 0, randX=0, randY=0;
  srand(time(NULL));
  for (size_t i = 0; i < nonZeros; i+=2) {
    do {
      randX = rand() % size;
      randY = rand() % size;
    } while(randX == randY);
    *(matrix+randX*size+randY) = 1;
    *(matrix+randY*size+randX) = 1;
    sum += 2;
  }
  double sparsity = 1.0 - ((double)sum)/((double)(size*size));

  int save = SAVE;
  if ( save){
  FILE* file = fopen("input.txt", "w");
      if(file == NULL)
        exit(0);

      for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < size; j++) {
          fprintf(file, "%d, ", *(matrix+i*size+j));
        }
        fprintf(file, "\n");
      }
      fclose(file);
  }
}


void ReverseCuthilMckee(int *matrix , int *degreeArray ,int *  R){
unsigned cap = QUEUE_SIZE;
Queue *Q  = createQueue(cap);
int Rsize = 0;

int *notAdded = (int*)malloc(N*sizeof(int));
if( notAdded == NULL){
  printf("Error at memory Initialization \n");
  exit(0);
}
int currentIndex ;
omp_lock_t writelock;
omp_init_lock(&writelock);

int i ;
#pragma omp parallel for  schedule(static,CHUNKSIZE) num_threads (MAX_THREADS)  private(i) shared(notAdded)
  for ( i = 0; i < N; i++) {
    notAdded[i] = 1;
  }

//Iterate for all the rows
while(Rsize != N){

  int minDegreeINdex  = findMinIdxParallel(degreeArray, notAdded , N, MAX_THREADS);


  enqueue( Q, minDegreeINdex);
  notAdded[minDegreeINdex] = 0 ;

  while(!(isEmpty(Q))){

   currentIndex = dequeue(Q);
   //Find the neighbors of the current Index
   int *neighbors = (int *)malloc(degreeArray[currentIndex]*sizeof(int));
   if(neighbors == NULL){
     printf("Error at memory Initialization \n");
     exit(0);
   }
   int neighborsCounter = 0;

#pragma omp parallel for  schedule(static,CHUNKSIZE) num_threads (MAX_THREADS)  private(i) shared(notAdded , matrix)
   for (int i = 0; i < N ; i++){
     if( i!= currentIndex && *(matrix+(currentIndex)*N + i) == 1 && notAdded[i] == 1 ){
       omp_set_lock(&writelock);
       neighbors[neighborsCounter++] = i;
       notAdded[i] = 0;
       omp_unset_lock(&writelock);
     }
   }

  //sort the neighbors found by their degree value in order to get the increasing order
  int *tempArray = (int *)malloc(neighborsCounter*sizeof(int));

//Calculate a temp Array that has the indexes of the neighbor nodes
//In order to sort the degreeArray (increasing order)
  for(int i = 0; i < neighborsCounter ; i++){
    tempArray[i] = degreeArray[neighbors[i]];
  }

  mergeSort(degreeArray , tempArray  , 0, neighborsCounter - 1  );


  #pragma omp parallel for ordered shared(Q,neighborsCounter) private(i)
           for( i = 0; i < neighborsCounter ; i++){
             #pragma omp ordered
             {
                    enqueue( Q, neighbors[i]);
             }
           }


  R[Rsize++] = currentIndex;

  free(neighbors);
  }
 }
  free(Q);

  //Reversing
  int nSize  = N;

  if (nSize % 2 == 0)
      nSize -= 1;

     nSize= nSize / 2;

      for (int i = 0; i <= nSize; i++) {
          int j = R[N - 1 - i];
          swapElement(&R[N - 1 - i] , &R[i]);
          R[i] = j;
        }

  }

  int* matrixReorder(int* matrix, int* R, int size) {
    int* newMatrix =(int*) calloc(size*size, sizeof(int));
    if(newMatrix == NULL){
      printf("Error at memory Initialization \n");
      exit(0);
    }

    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        if(*(matrix+size*R[i]+j) == 1) {
          for (int k = 0; k < size; k++) {
            if(R[k] == j) {
              *(newMatrix+size*i+k) = 1;
            }
          }
        }
      }
    }

    return newMatrix;
  }

void sortByDegree(int *neighbors ,int  *degreeArray , int neighborsCounter){

int *tempArray = (int *)malloc(neighborsCounter*sizeof(int));
if(tempArray == NULL){
  printf("Error at memory Initialization \n");
  exit(0);
}

for(int i = 0; i < neighborsCounter ; i++){
  tempArray[i] = degreeArray[neighbors[i]];
}
//Sort the array of the neighbors in increasing order by their degree
for (int i = 0; i < neighborsCounter-1; i++){
    for (int j = 0; j < neighborsCounter-i-1; j++){
      if (tempArray[j] > tempArray[j+1]){
          swapElement(&tempArray[j], &tempArray[j+1]);
          swapElement(&neighbors[j], &neighbors[j+1]);
      }
    }
  }
}


void swapElement(int* one, int* two) {
  int temp = *one;
  *one = *two;
  *two = temp;
}
  int minDegreeINdex = 0 ;

  void degreeCalculator(int* matrix, int * degreeArray) {

    int i;

    #pragma omp parallel num_threads(MAX_THREADS) private(i)
    {
      int sum=0;
      #pragma omp for schedule(dynamic , CHUNKSIZE)
        for ( i = 0; i < N; i++) {
          for (size_t j = 0; j < N; j++) {
            sum += *(matrix+ i*N+j);
          }
          degreeArray[i] = sum;
          sum = 0;
        }
    }

  }

int findMinIdxParallel(int* degrees, int* notVisited, int size, int threads_num) {
  int degArray[threads_num];
  int idxArray[threads_num];
  int tid;

  // find the min degree
  int i;
  int CHUNKSIZE = size / threads_num;
  #pragma omp parallel num_threads(MAX_THREADS) private(i, tid)
  {
    tid = omp_get_thread_num();
    int minDegreeIndex = -1;  // The pos of min degree node in matrix
    int minDegree = size+10; // A node can not have degree > SIZE
    #pragma omp for schedule(static, CHUNKSIZE)
      for (i = 0; i < size; i++) {
        if(degrees[i] < minDegree && notVisited[i] == 1) {
          minDegreeIndex = i;
          minDegree = degrees[i];
        }
      }
      degArray[tid] = minDegree;
      idxArray[tid] = minDegreeIndex;
  }

  int globalMinIdx = idxArray[0];
  int globalMinDeg = degArray[0];

  for (size_t i = 1; i < threads_num; i++) {
    if(degArray[i] < globalMinDeg) {
      globalMinDeg = degArray[i];
      globalMinIdx = idxArray[i];
    }
  }

  return globalMinIdx;
}

void merge(int * arr,int * idx, int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    //Temp arrays for the Left and Right sub problem
    int L[n1], R[n2];
    int Lidx[n1], Ridx[n2];

    //Copying the data to the subarrays using two different threads , one for each
    #pragma omp parallel sections  shared(arr,idx)
    {
      #pragma omp section
      for (i = 0; i < n1; i++){
          L[i] = arr[l + i];
          Lidx[i] = idx[l + i];
      }

      #pragma omp section
      for (j = 0; j < n2; j++){
          R[j] = arr[m + 1 + j];
          Ridx[j] = idx[m + 1 + j];
      }
    }
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k]  = L[i];
            idx[k] = Lidx[i];
            i++;
        }
        else {
            arr[k] = R[j];
            idx[k] = Ridx[j];
            j++;
        }
        k++;
    }

    while (i < n1) {
        arr[k] = L[i];
        idx[k] = Lidx[i];
        i++;
        k++;
    }

    while (j < n2) {
        arr[k] = R[j];
        idx[k] = Ridx[j];
        j++;
        k++;
    }
}

// l is for left index and r is right index of the
   sub-array of arr to be sorted //
void mergeSort(int * arr,int * idx, int l, int r)
{
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        int m = l + (r - l) / 2;


          // one thread sorts the first half, and another the second
          #pragma omp parallel sections  shared(arr,idx)
          {
            #pragma omp section
            mergeSort(arr,idx, l, m);
            #pragma omp section
            mergeSort(arr,idx, m + 1, r);
          }

            merge(arr,idx, l, m, r);


    }
}


void outputWrite(int* matrix, int rows, int col, char* file_path)
{
  FILE* file = fopen(file_path, "w");
  if(file == NULL)
    exit(0);

  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < col; j++) {
      fprintf(file, "%d, ", *(matrix+i*col+j));
    }
    fprintf(file, "\n");
  }
  fclose(file);
}
