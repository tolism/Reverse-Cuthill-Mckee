/*
* FILE: sequential.c
* THMMY, 7th semester, Parallel and Distributed Systems: 4th assignment
* Sequential Implementation of the Reverse Cuthill Mckee Algorithm
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



#define N 1000
#define SPARSITY_LIMIT 0.85
#define QUEUE_SIZE 10000
#define MODE 2
#define SAVE 0



//Function Decleration
void MatrixInitialization(int * matrix , int size , double sparsityLimit , int mode );
void degreeCalculator(int *matrix , int * degreeArray);
void ReverseCuthilMckee(int *matrix , int *degreeArray , int * R );
void swapElement(int *one, int  *two);
void sortByDegree(int *neighbors ,int  *degreeArray , int neighborsCounter);
int* matrixReorder(int* matrix, int* R, int size);
void outputWrite(int* matrix, int rows, int col, char* file_path);


int main(int  argc, const char* argv[]){


int *  matrix = (int * )calloc(N*N,sizeof(int));
int * degreeArray = (int *)calloc(N,sizeof(int));
int  * R = (int *)calloc(N,sizeof(int));

if(matrix == NULL || degreeArray == NULL || R == NULL){
  printf("Error at memory Initialization \n");
  exit(0);
}

struct timeval start, end;

MatrixInitialization(matrix, N , SPARSITY_LIMIT , MODE);



gettimeofday(&start, NULL);
degreeCalculator(matrix , degreeArray) ;
ReverseCuthilMckee(matrix , degreeArray , R);
gettimeofday(&end, NULL);
int numberOfN = N;
double time = ((double)((end.tv_sec*1e6 + end.tv_usec) - (start.tv_sec*1e6 + start.tv_usec)))*1e-6;
printf("Execution Time: %lf sec sequential version N=%d\n", time , numberOfN);

int save = SAVE;
if ( save){
   int* reorderedMatrix = matrixReorder( matrix, R, N) ;
   outputWrite(reorderedMatrix, N, N, "outputMatrix.txt");
}
};


//Matrix Initialization function
void MatrixInitialization (int * matrix , int size , double sparsityLimit , int mode  ){


  if(mode == 1) {
      FILE* file = fopen("input.txt", "r");
      if(file == NULL)
        exit(0);

      int i=0, j=0;
      int value;
      char* ch = (char*) malloc(sizeof(char));


      while(1) {
        *ch = fgetc(file);
        if(*ch == '0' || *ch == '1') {
          value = atoi(ch);
          *(matrix+i*size+j) = value;
          j++;
          if(j >= size){
            i++;
            j=0;
          }
        }

        if(*ch == EOF) {
          break;
        }

      }
      fclose(file);
    }
else{

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



  int save = SAVE ;
  if (save){
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
}

//THe Reverse Cuthil Mckee function
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


  for (int  i = 0; i < N; i++) {
    notAdded[i] = 1;
  }

//Iterate for all the rows
while(Rsize != N){

  int minDegreeINdex = 0 ;
  for (int i = 0; i< N ; i++){
    if(degreeArray[i] < degreeArray[minDegreeINdex] && notAdded[i]==1){
      minDegreeINdex = i;
    }
  }

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


   for (int i = 0; i < N ; i++){
     if( i!= currentIndex && *(matrix+(currentIndex)*N + i) == 1 && notAdded[i] == 1 ){
       neighbors[neighborsCounter++] = i;
       notAdded[i] = 0;
     }
   }

  //sort the neighbors found by their degree value in order to get the increasing order
  sortByDegree(neighbors , degreeArray , neighborsCounter);

  for(int i = 0; i < neighborsCounter ; i++){
      enqueue( Q, neighbors[i]);
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

//Sort the neighbors by their degree
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

//Swwapping element function
void swapElement(int* one, int* two) {
  int temp = *one;
  *one = *two;
  *two = temp;
}

//Implementation of the Degree Calculator function
void degreeCalculator(int *matrix,int *  degreeArray){
  int sum = 0 ;
  for(int i = 0; i< N; i++){
    for(int j = 0; j < N ; j++){
      sum +=*(matrix + i*N + j);
    }
    degreeArray[i]=sum;
    sum = 0;
  }
}

//Reordering the matrix using R
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

//Function used to write the reversed matrix to a file
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
