//------------------------------------------------------------------------- 
// This is supporting software for CS415/515 Parallel Programming.
// Copyright (c) Portland State University
//------------------------------------------------------------------------- 

// Jacobi method for solving a Laplace equation.  
//
// Usage: ./jacobi [N]
// 
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPSILON 0.001 	// convergence tolerance
#define VERBOSE 0 	// printing control

// Initialize the mesh with a fixed set of boundary conditions.
// 
void init_array(int n, double a[n][n])  {
  int i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) 
      a[i][j] = 0;
  }
  for (i = 1; i < n-1; i++) {
    a[n-1][i] = 1.0;
    a[i][n-1] = 1.0;
  }
}

// Display the whole mesh.
// 
void print_array(int n, double a[n][n])  {
  int i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      printf("%8.4f ", a[i][j]);
    printf("\n");
  }
}

// Jacobi iteration -- return the iteration count.
// 
int jacobi(int n, double x[n][n], double epsilon) {
  double xnew[n][n];	// buffer for new values    
  double delta;		// measure of convergence   
  int cnt = 0;		// iteration counter              
  int i, j;

  do {	
    delta = 0.0;
    for (i = 1; i < n-1; i++) {
      for (j = 1; j < n-1; j++) {
	xnew[i][j] = (x[i-1][j] + x[i][j-1] + x[i+1][j] + x[i][j+1]) / 4.0;
	delta = fmax(delta, fabs(xnew[i][j] - x[i][j]));
      }
    }	
    for (i = 1; i < n-1; i++) {
      for (j = 1; j < n-1; j++) {
	x[i][j] = xnew[i][j];
      }
    }
    cnt++;
    if (VERBOSE) {
      printf("Iter %d: (delta=%6.4f)\n", cnt, delta);
      print_array(n, x);
    }
  } while (delta > epsilon);
  return cnt;
}

//gauss_seidel version
int gauss_seidel(int n, double x[n][n], double epsilon) {
  double delta, temp;
  int cnt = 0;
  int i, j;

  do {
    delta = 0.0;
    for (i = 1; i < n-1; i++) {
      for (j = 1; j < n-1; j++) {
        temp = x[i][j];
        x[i][j] = (x[i-1][j] + x[i][j-1] + x[i+1][j] + x[i][j+1]) / 4.0;
	delta = fmax(delta, fabs(x[i][j] - temp));
      }
    }
    cnt++;
    if (VERBOSE) {
      printf("Iter %d: (delta=%6.4f)\n", cnt, delta);
      print_array(n, x);
    }
  } while (delta > epsilon);
  return cnt;
}

int red_black(int n, double x[n][n], double epsilon) {
  double delta, delta_r, delta_b, temp_r, temp_b;
  int cnt = 0;
  int i, j;

  do {
    delta_r = 0.0;
    for (i = 1; i < n-1; i++) {
      for (j = 1; j < n-1; j++) {
        if (!((i + j) % 2)){
          temp_r = x[i][j];
          x[i][j] = (x[i-1][j] + x[i][j-1] + x[i+1][j] + x[i][j+1]) / 4.0;
	  delta_r = fmax(delta_r, fabs(x[i][j] - temp_r));
        }
      }
    }

    delta_b = 0.0;    
    for (i = 1; i < n-1; i++) {
      for (j = 1; j < n-1; j++) {
        if ((i + j) % 2){
          temp_b = x[i][j];
          x[i][j] = (x[i-1][j] + x[i][j-1] + x[i+1][j] + x[i][j+1]) / 4.0;
	  delta_b= fmax(delta_b, fabs(x[i][j] - temp_b));
        }
      }
    }

    delta = fmax(delta_b, delta_r);
    cnt++;

    if (VERBOSE) {
      printf("Iter %d: (delta=%6.4f)\n", cnt, delta);
      print_array(n, x);
    }
  } while (delta > epsilon);
  return cnt;
}

// Main routine.
//
int main(int argc, char **argv) {

  int n = 32;  	   	// mesh size, default 8 x 8
  if (argc > 1) {  	// check command line for overwrite
    if ((n = atoi(argv[1])) < 2) {
      printf("Mesh size must must be greater than 2, use default\n");
      n = 8;
    }
  }

  double a[n][n];	// mesh array
  init_array(n, a);

  // Jacobi iteration, return value is the total iteration number
  int j_cnt = jacobi(n, a, EPSILON);
  printf("Jacobi:\n");
  printf("Mesh size: %d x %d, epsilon=%6.4f, total Jacobi iterations: %d\n", 
	 n, n, EPSILON, j_cnt);
  if (VERBOSE) 
    print_array(n, a);

//  print_array(n, a);

  // reinitialize the array
  init_array(n, a);
  // Gauss-Seidel method, return is iteration count
  int gs_cnt = gauss_seidel(n, a, EPSILON);
  printf("Gauss-Seidel:\n");
  printf("Mesh size: %d x %d, epsilon: %6.4f, iterations: %d\n", 
          n, n, EPSILON, gs_cnt);

//  print_array(n, a);

  //reinitialize the array
  init_array(n, a);
  // Gauss-Seidel red/black method, return is iteration count
  int rb_cnt = red_black(n, a, EPSILON);
  printf("Red/Black:\n");
  printf("Mesh size: %d x %d, epsilon: %6.4f, iterations: %d\n",
          n, n, EPSILON, rb_cnt);

//  print_array(n, a);
}
