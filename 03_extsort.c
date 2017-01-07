//------------------------------------------------------------------------- 
// Copyright (c) Thomas Van Klaveren 2015 
//------------------------------------------------------------------------- 

//  This program takes a byte file of N unsorted integers, supplied
//  by the user at runtime and 
//  sorts the integers using mpi commands. The result is output to a 
//  file spiceified by the user at runtime.  
//  Application assumes that N is greater than 10P where P = number of 
//  processes.
//
// Usage: 
//   linux> mpirun -hostflie <hostfile> -n <#processes> extsort 
//          [<inputfile> <outputfile>]
// 
// 
#define _BSD_SOURCE
#include <unistd.h>	// for gethostname()
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#define TAG 1001
#define MINSIZE 10

// find min value
double min(double *array, int length){
  double s = array[0];

  for (int i = 1; i < length; i++){
    if (array[i] < s){
      s = array[i];
    }
  }
 
  return s;
}

//find max value
double max(double *array, int length){
  double m = array[0];

  for (int i = 1; i < length; i++){
    if (array[i] > m){
      m = array[i];
    }
  }

  return m;
}
  

// Swap two array elements 
//
void swap(int *array, int i, int j) {
  if (i == j) return;
  int tmp = array[i];
  array[i] = array[j];
  array[j] = tmp;
}

// Bubble sort for the base cases
//
void bubblesort(int *array, int low, int high) {
  if (low >= high) 
    return;
  for (int i = low; i <= high; i++)
    for (int j = i+1; j <= high; j++) 
      if (array[i] > array[j])
	swap(array, i, j);
}

// Pick an arbitrary element as pivot. Rearrange array 
// elements into [smaller one, pivot, larger ones].
// Return pivot's index.
//
int partition(int *array, int low, int high) {
  int pivot = array[high]; 	// use highest element as pivot 
  int middle = low;
  for(int i = low; i < high; i++)
    if(array[i] < pivot) {
      swap(array, i, middle);
      middle++;
    }
  swap(array, high, middle);
  return middle;
}
 
// QuickSort an array range
// 
void quicksort(int *array, int low, int high) {
  if (high - low < MINSIZE) {
    bubblesort(array, low, high);
    return;
  }
  int middle = partition(array, low, high);
  if (low < middle)
    quicksort(array, low, middle-1);
  if (middle < high)
    quicksort(array, middle+1, high);
}
 

int main(int argc, char *argv[])
{
  int nprocs, rank;
  char host[50];
  MPI_Status status;
  MPI_Offset filesize;
  int N, my_count;
  MPI_File in, out;
  int *array, *pivot, *count, *my_bucket, *bucket_count;
  int **bucket;
  double begin, end, end_no_io;
  double begin_read, end_read, begin_write, end_write;
  double *time_start, *time_io, *time_no_io;

  //get time including io
  begin = MPI_Wtime();
//  printf("beg w/ io at %f\n", begin);

  if (argc != 3) {
    printf("Useage: ./file <input> <output>\n");
    exit(1);
  }

  gethostname(host, 50);

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if (nprocs < 2) {
    printf("Need at least 2 processes.\n");
    MPI_Finalize();
    return(1);
  }

//  printf("P%d/%d started on %s ...\n", rank, nprocs, host);

  begin_read = MPI_Wtime();

  if (rank == 0) {

    //get size of file and number of elements,
    //read data from input file to array
    MPI_File_open(MPI_COMM_SELF, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
    MPI_File_get_size(in, &filesize);
    N = (filesize/sizeof(int));
    array = (int *)(malloc(sizeof(int) * N));
    MPI_File_read(in, array, N, MPI_INT, &status);
    MPI_File_close(&in);
  }

  end_read = MPI_Wtime();

  if (rank == 0) {

    // set up time array
    time_start = (double *)(malloc(sizeof(double) * nprocs));
    time_io = (double *)(malloc(sizeof(double) * nprocs));
    time_no_io = (double *)(malloc(sizeof(double) * nprocs));

    //sort first 10P elements
    quicksort(array, 0, (10 * nprocs)-1);

    //store elements at position 10, 20, 30... (10*(P-1)) as pivot vals
    pivot = (int *)(malloc(sizeof(int) * nprocs));
    for (int i = 0; i < (nprocs -1); i ++){
      pivot[i] = array[(i + 1) * 10];
    }
    pivot[nprocs - 1] = N;

    //get counts for bucket sizes
    count = (int *)(malloc(sizeof(int) * nprocs));
    count[0] = pivot[0];
    for (int i = 1; i < nprocs; i++){
      count[i] = pivot[i] - pivot[i-1];
    }

  }

  //parallel computation starts here
  //scatter count so everyone knows their bucket size
  MPI_Scatter(count, 1, MPI_INT, &my_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
//  printf("node %d/%d got count = %d\n", rank, nprocs, my_count);

  //give each proc their bucket
  my_bucket = (int *)(malloc(sizeof(int) * my_count));

  if (rank == 0){
    //put data into appropriate bucket
    bucket_count = (int *)(malloc(sizeof(int) * nprocs));
    for (int i = 0; i < nprocs; i++){
      bucket_count[i] = 0;
    }

    bucket = (int **)(malloc(sizeof(int *) * nprocs));
    for (int i = 0; i < nprocs; i++){
      bucket[i] = (int *)(malloc(sizeof(int) * count[i]));
    }

    for (int i = 0; i < N; i++){
      for (int j = 0; j < nprocs; j++){
        if (array[i] <= pivot[j]){
          bucket[j][bucket_count[j]] = array[i];
          bucket_count[j]++;
          break;
        }
      }
    }

    //free array, pivot, bucket_count memory
    free(array);
    free(pivot);
    free(bucket_count);

    //send buckets to other nodes
    for (int i = 1; i < nprocs; i++){
      MPI_Send(bucket[i], count[i], MPI_INT, i, TAG, MPI_COMM_WORLD);
    }

    //note: this is wasteful - could just qsort on bucket[root]
    //however, decided to linearly move to my_bucket so that
    //the concurrent code below would encompass root proc.
    for(int i = 0; i < my_count; i++){
      my_bucket[i] = bucket[rank][i];
    }

    //free bucket and count memory
    for (int i = 0; i < nprocs; i++){
      free(bucket[i]);
    }

    free(bucket);
    free(count);

  }

  else{
    //give the other procs their buckets
    MPI_Recv(my_bucket, my_count, MPI_INT, 0, TAG, MPI_COMM_WORLD, &status);
//    printf("node %d/%d got a bucket\n", rank, nprocs);
  }

  //parallel computation going forward
  //qsort the bucket
  quicksort(my_bucket, 0, my_count-1);
  begin_write = MPI_Wtime();

  //write my section to the file
  //open file
  MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_CREATE|MPI_MODE_RDWR, 
                MPI_INFO_NULL, &out);

  //go to the appropriate space in the file
  MPI_File_set_view(out, ((my_bucket[0] -1) * sizeof(int)), MPI_INT, 
                    MPI_INT, "native", MPI_INFO_NULL);

  //write in my section of the file and close file
  MPI_File_write(out, my_bucket, my_count, MPI_INT, &status);

  MPI_File_close(&out);

  end_write = MPI_Wtime();

  //free my_bucket memory
  free(my_bucket);

  //get time for each proc taking io into account
  end = MPI_Wtime();
//  printf("end w/io at %f no io at %f\n", end, end - (end_read - begin_read) -
//         (end_write - begin_write));

  //organize and deliver time to root
  end_no_io = end - (end_read - begin_read) - (end_write - begin_write);

  MPI_Gather(&end, 1, MPI_DOUBLE, 
             time_io, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Gather(&end_no_io, 1, MPI_DOUBLE, 
             time_no_io, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Gather(&begin, 1, MPI_DOUBLE, 
             time_start, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank == 0){
    //calculate final timing data for testing
    begin = min(time_start, nprocs);
    end = max(time_io, nprocs);
    end_no_io = max(time_no_io, nprocs);
  printf("begin = %f, end = %f, end no io = %f\n", begin, end, end_no_io);
    printf("time including IO: %f\n"
           "time excluding IO: %f\n",
           end - begin, end_no_io - begin);

    free(time_start);
    free(time_io);
    free(time_no_io);
  }
  MPI_Finalize();


}
