#ifndef QUEUE_H
#define QUEUE_H


typedef struct  {
    int first, last, size;
    unsigned capacity;
    int* array;
}Queue;


Queue* createQueue(unsigned capacity);
int isFull( Queue* queue);
int isEmpty( Queue* queue);
void enqueue( Queue* queue, int item);
int dequeue(Queue* queue);
int first(Queue* queue);
int last( Queue* queue);

#endif
