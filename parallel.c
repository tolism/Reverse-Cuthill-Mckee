#include <stdio.h>
#include <stdlib.h>
#include "queue.h"
#include <time.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <sys/time.h>




// !@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// NA VAL ELEGXO SE OLES TIS PIPES THS MALLOC
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define N 20000
#define SPARSITY_LIMIT 0.7
#define QUEUE_SIZE 100000
#define MAX_THREADS 4

void MatrixInitialization(int * matrix , int size , double sparsityLimit);
void degreeCalculator(int *matrix , int * degreeArray);
void ReverseCuthilMckee(int *matrix , int *degreeArray , int * R );
void swapElement(int *one, int  *two);
void sortByDegree(int *neighbors ,int  *degreeArray , int neighborsCounter);
void matrixReorder(int * matrix , int * R );



int main(int  argc, const char* argv[]){

int *  matrix = (int * )malloc(sizeof(int)*N*N);
struct timeval start, end;

MatrixInitialization(matrix, N , SPARSITY_LIMIT);
//
//
// for(int i = 0 ; i < N ; i++){
//   for(int j = 0 ; j < N ; j++){
//     printf("%d" , *(matrix + i*N + j));
//   }
//   printf("\n");
// }


int * degreeArray = (int *)malloc(N*sizeof(int));
degreeCalculator(matrix , degreeArray) ;


int  * R = (int *)malloc(N*sizeof(int));

gettimeofday(&start, NULL);
ReverseCuthilMckee(matrix , degreeArray , R);
gettimeofday(&end, NULL);

double time = ((double)((end.tv_sec*1e6 + end.tv_usec) - (start.tv_sec*1e6 + start.tv_usec)))*1e-6;
printf(" >>> Execution Time: %lf sec\n", time);

//Printing the results
// for(int i = 0 ; i < N; i++){
//   printf("%d " , R[i]);
// }
// printf("\n");

// for(int i = 0 ; i < N ; i++){
//   for(int j = 0 ; j < N ; j++){
//     printf("%d" , *(matrix + R[i]*N + j));
//   }
//   printf("\n");
// }


};

//Matrix Initialization function
void MatrixInitialization (int * matrix , int size , double sparsityLimit  ){

  //
  // FILE* file = fopen("input.txt", "r");
  //     if(file == NULL)
  //       exit(0);
  //
  //     int i=0, j=0;
  //     int value;Cilk_lockvarmylock;
  //     char ch;
  //
  //
  //     while(1) {
  //       ch = fgetc(file);
  //       value = atoi(&ch);
  //       if(ch == '0' || ch == '1') {
  //         *(matrix+i*N+j) = value;
  //         j++;
  //         if(j >= N){
  //           i++;
  //           j=0;
  //         }
  //       }
  //
  //       if(ch == EOF) {
  //         break;
  //       }
  //
  //     }
  //     fclose(file);

  //srand(time(NULL));
    double randomNumber=0;

    while(1){
      double sum=0;
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          randomNumber = (rand() % 100);
        //  printf("%f  " , randomNumber);
          if(randomNumber > 70 ){
            randomNumber = 1;
            sum += randomNumber;
          }
          else{
            randomNumber = 0;
          }
          *(matrix+i*N+j) = (int)randomNumber;
          *(matrix+j*N+i) = (int)randomNumber;


        }
      }
      double sparsity = (((double)(N*N)) - sum) / (double)(N*N);
  //    printf(" >> sparsity: %lf\n", sparsity);
      if(sparsity > SPARSITY_LIMIT ){
        break;
      }
    }
}

void ReverseCuthilMckee(int *matrix , int *degreeArray ,int *  R){
unsigned cap = QUEUE_SIZE;
Queue *Q  = createQueue(cap);
int Rsize = 0;
int *notAdded = (int*)malloc(N*sizeof(int));
int currentIndex ;


  cilk_for (int  i = 0; i < N; i++) {
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
   int neighborsCounter = 0;


   cilk_for (int i = 0; i < N ; i++){
     if( i!= currentIndex && *(matrix+(currentIndex)*N + i) == 1 && notAdded[i] == 1 ){
       neighbors[neighborsCounter++] = i;
       notAdded[i] = 0;
     }
   }

  //sort the neighbors found by their degree value in order to get the increasing order
  sortByDegree(neighbors , degreeArray , neighborsCounter);

  cilk_for(int i = 0; i < neighborsCounter ; i++){
      enqueue( Q, neighbors[i]);
  }
  R[Rsize++] = currentIndex;

  free(neighbors);
  }
 }
  free(Q);

  //To check if reverse works fine
  // for(int i = 0 ; i < N; i++){
  //   printf("%d " , R[i]);
  // }
  // printf("\n");
  // }


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

void matrixReorder(int * matrix , int * R ){











}


void sortByDegree(int *neighbors ,int  *degreeArray , int neighborsCounter){

int *tempArray = (int *)malloc(neighborsCounter*sizeof(int));

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
