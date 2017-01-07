//Thomas Van Klaveren assignment 1
//Parallel Programming
//A global array of size N contains the integers to be sorted
//A global task queue is initialized with the sort range [0, N-1]


#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sched.h>
#include <unistd.h>
#include <time.h>

//threshhold for switching to bubble sort
#define MINSIZE 10 



//----------------------------------------------------------------------
//Task and queue representaions
typedef struct task_ {
  int low;
  int high;
  struct task_ *next;
} task_t;

typedef struct queue_{
  task_t *head;
  task_t *tail;
  int length;
} queue_t;

//create a new task
task_t *create_task(int low, int high){
  task_t *task = (task_t *) malloc(sizeof(task_t));
  task->low = low;
  task->high = high;
  task->next = NULL;
  return task;
}

//initialize the task queue
queue_t *init_queue(){
  queue_t *queue = (queue_t *) malloc(sizeof(queue_t));
  queue->head = queue->tail = NULL;
  queue->length = 0;
  return queue;
}

//add a task to tail of queue
void add_task(queue_t *queue, task_t *task){
  if(!queue->tail){
    queue->head = queue->tail = task;
  }
  else{
    queue->tail->next = task;
    queue->tail = task;
  }

  queue->length++;
}

//remove a task from the head of the queue (return NULL if empty)
task_t *remove_task(queue_t *queue){
  task_t *task = NULL;

  if(queue->length > 0){
    task = queue->head;

    if(queue->head == queue->tail){
      queue->head = queue->tail = NULL;
    }
    else{
      queue->head = queue->head->next;
    }

    queue->length--;
  }

  return task;
}
//-------------------------------------------------------------------------

//global shared variables (array, queue, count, mutex lock)
int *array = NULL;
queue_t *queue;
int count = 0;
int N = 0;
pthread_mutex_t queue_lock;
pthread_mutex_t sum_lock;
pthread_cond_t length_cond;

//print array for testing purposes
void print_array(int *array, int low, int high){
  printf("low = %d, a[%d] = %d\n"
         "high = %d, a[%d] = %d\n", 
          low+1, low+1, array[low], high, high, array[high-1]);
  for(int i = 0; i<N; i++){
    printf("%d ", array[i]);
  }
  printf("\n");
}

//swap 2 elements of an array
void swap(int *array, int i, int j){
  if (i == j){
    return;
  }

  int temp = array[i];
  array[i] = array[j];
  array[j] = temp;
}


//initialize an array of N elements
//1. geneate [1, 2, 3, ... N]
//2. perform a random permutation
int *init_array(int N){
  int *array = (int *) malloc(sizeof(int) * N);

  for (int i = 0; i < N; i++){
    array[i] = i + 1;
  }

  srand(time(NULL));
  for (int i = 0; i < N; i++){
    int j = (rand() * 1./RAND_MAX) * (N-1);
    swap(array, i, j);
  }

  printf("Initialized array to a random permutation of [1..%d]\n", N);
  return array;
}

//verify the result
void verify_array(int *array, int N){
  for (int i = 0; i < N-1; i++){
    if (array[i] > array[i+1]){
      printf("FAILED: array[%d] = %d, array[%d] = %d\n", 
             i, array[i], i+1, array[i+1]);
      return;
    }
  }

  printf("Result verified!\n");
}

//bubble sort for the base cases
void bubblesort(int *array, int low, int high){
  if(low >= high){
    return;
  }

  for(int i = low; i <= high; i++){
    for(int j = i+1; j <= high; j++){
      if(array[i] > array[j]){
        swap(array, i, j);
      }
    }
  }
}

//pick element as pivot
//rearrange array elements into [smaller elements, pivot, larger elements]
//retunr pivot's final index
int partition(int *array, int low, int high){
  int pivot = array[high];
  int middle = low;

  for(int i = low; i < high; i++){
    if(array[i] < pivot){
      swap(array, i, middle);
      middle++;
    }
  }

  swap(array, high, middle);
  return middle;
}


//quicksort an array
void quicksort(int *array, int low, int high){
  if(high - low < MINSIZE){
    bubblesort(array, low, high);

    //update global count 
    pthread_mutex_lock(&sum_lock);
    count += (high - low)+1;

    //add dummy task if global count > N
    if(count >= N){
      pthread_mutex_unlock(&sum_lock);
      task_t *dummy = NULL;
      dummy = create_task(0, 0);
      pthread_mutex_lock(&queue_lock);
      add_task(queue, dummy);

      //signal waiting condition since queue length increase
      pthread_cond_signal(&length_cond);
      pthread_mutex_unlock(&queue_lock);
    }

    else{
      pthread_mutex_unlock(&sum_lock);
    }

    return;
  }

  //partition the array
  int middle = partition(array, low, high);

  //update global count and add dummy task if count > N
  pthread_mutex_lock(&sum_lock);
  count++;

  if(count >= N){
    pthread_mutex_unlock(&sum_lock);
    task_t *dummy = NULL;
    dummy = create_task(0, 0);
    pthread_mutex_lock(&queue_lock);
    add_task(queue, dummy);

    //signal waiting condition since queue length increase
    pthread_cond_signal(&length_cond);
    pthread_mutex_unlock(&queue_lock);
  }

  else{
    pthread_mutex_unlock(&sum_lock);
  }

  if (low < middle){
    //create task and add to queue for next avail thread
    //for the array elements on the left side of the partition
    task_t *task = NULL;
    task = create_task(low, middle-1);
    pthread_mutex_lock(&queue_lock);
    add_task(queue, task);
    //signal waiting threads to wake up for new task on queue
    pthread_cond_signal(&length_cond);
    pthread_mutex_unlock(&queue_lock);
  }


  if (middle < high){
    //recursively quicksort on the elements 
    //on the right side of the partition
    quicksort(array, middle+1, high); 
  }
}


//worker routine that each thread will call
void worker(long wid){
  printf("worker %ld started on %d\n", wid, sched_getcpu());
  int l_count = 0;
  task_t *task;

  //the following will repeat if the count is less than
  //the array size at the end of the block
  do{
    pthread_mutex_lock(&sum_lock);
    l_count = count;
    pthread_mutex_unlock(&sum_lock);

    //if the local count is less than the array size
    //access the queue and wait for a task if the 
    //queue is empty
    if(l_count < N){
      pthread_mutex_lock(&queue_lock);
      while(queue->length < 1){
        pthread_cond_wait(&length_cond, &queue_lock); 
      }

      //after wake up from waiting, get task and quicksort
      task = remove_task(queue);
      pthread_mutex_unlock(&queue_lock);
      quicksort(array, task->low, task->high);
    }

    pthread_mutex_lock(&sum_lock);
    l_count = count;
    pthread_mutex_unlock(&sum_lock);
  }while (l_count < N);
     
}



//main routine
int main(int argc, char **argv){


  //default number of threads
  int num_thread = 1;

  //check user inputs
  if(argc < 2){
    printf("Usage:  ./qsortpthrd <N> <num_thread>\n");
    exit(0);
  }

  if((N = atoi(argv[1])) < 2){
    printf("<N> must be greater than 2\n");
    exit(0);
  }

  if((num_thread = atoi(argv[2])) < 1){
    printf("<num_thread> must be greater than 0\n");
    exit(0);
  }


  //initialize array, queue, locks
  array = init_array(N);

  queue = init_queue();

  pthread_t thread[num_thread];

  pthread_mutex_init(&queue_lock, NULL);
  pthread_mutex_init(&sum_lock, NULL);
  pthread_cond_init(&length_cond, NULL);

  //create first task
  task_t *task_one = NULL;
  task_one = create_task(0, N-1);
  pthread_mutex_lock(&queue_lock);
  add_task(queue, task_one);
  pthread_mutex_unlock(&queue_lock);


  //create numThreads-1 worker threads to execute worker()
  for (long k = 0; k < num_thread - 1; k++){
    pthread_create(&thread[k], NULL, (void*)worker, (void*)k);
  }
  
  //make main thread run worker()
  worker(num_thread - 1);
  printf("%d threads created\n", num_thread);

  //main thread waits for worker threads to join
  for (long k = 0; k < num_thread - 1; k++){
    pthread_join(thread[k], NULL);
  }

  //varify the result
  verify_array(array, N);

return 0;
}
