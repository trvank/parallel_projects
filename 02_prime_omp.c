//------------------------------------------------------------------------- 
// This is supporting software for CS415/515 Parallel Programming.
// Copyright (c) Portland State University
//------------------------------------------------------------------------- 

// A sequential prime-finding algorithm.
//
// Usage: ./prime <N>
// 
//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

int main(int argc, char **argv) {
  long N;
  int num_thread;

  /* check command line first */
  if (argc < 3) {
    printf ("Usage: ./prime_omp <N> <num_thread>\n");
    exit(0);
  }
  if ((N=atoi(argv[1])) < 2) {
    printf ("N must be greater than 1\n");
    exit(0);
  }
  if ((num_thread = atoi(argv[2])) < 1){
    printf("<num_thread> must be greater than 0\n");
    exit(0);
  }

  omp_set_num_threads(num_thread);

#ifdef DEBUG
  printf("Finding primes in range 1..%d\n", N);
#endif

  long *array = (long *) malloc(sizeof(long) * (N+1));

  #pragma omp parallel for
  for (int i = 2; i <= N; i++){
    array[i] = 1;
//    printf("updated array on thread %d\n", omp_get_thread_num());
  }

  int limit = (int) sqrt((double) N);

//  #pragma omp parallel
  

//  #pragma omp for
//  #pragma omp single
  for (int i = 2; i <= limit; i++) {
//    #pragma omp task firstprivate(i)
    if (array[i] == 1) {
//      printf("cancelling mults of %d on thread %d\n", i, omp_get_thread_num());
      #pragma omp parallel for firstprivate(i)
//      #pragma omp single firsteprivate(i)
      for (int j = i+i; j <= N; j += i){
//        #pragma omp task firstprivate(j) 
	array[j] = 0;
//        printf("found nonprime of %d on thread %d\n", i, omp_get_thread_num());
      }
    }
  }

  

  int cnt = 0;
  omp_lock_t cnt_lock;
  omp_init_lock(&cnt_lock);
/*for(int i = 2; i < N+1; i++){
  printf("%d ", array[i]);
}
printf("\n\n");
*/
  #pragma omp parallel for reduction(+:cnt)
  for (int i = 2; i <= N; i++) {
    if (array[i] == 1){
//      omp_set_lock(&cnt_lock);
//      #pragma omp flush(cnt)
      cnt++;
//      printf("found prime on thread %d\n", omp_get_thread_num());
//      #pragma omp flush(cnt)
//      omp_unset_lock(&cnt_lock);
    }
  }

  printf("Total %d primes found\n", cnt);
}

