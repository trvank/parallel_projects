//------------------------------------------------------------------------- 
// This is supporting software for CS415/515 Parallel Programming.
// Copyright (c) Portland State University
//------------------------------------------------------------------------- 

// Jacobi method for solving a Laplace equation.  
//
// Usage: ./jacobi-shm -nl <#locales>
// 
//

use BlockDist;

config const epsilon = 0.001;	// convergence tolerance
config const verbose = false; 	// printing control
config const n = 8; 	        // mesh size (including boundary)

// Jacobi iteration -- return the iteration count.
// 
proc jacobi(D: domain(2), x: [D] real, epsilon: real) { 
  const ID = D.expand(-1,-1); 	// domain for interior points
  var xnew: [D] real;           // buffer for new values
  var delta: real; 		// measure of convergence 
  var cnt = 0;			// iteration counter

  var xcheck: [D] real = 8.8;   // to check if mapped correctly
  const BD = ID dmapped Block(boundingBox=ID); //mapping

  do {
    forall ij in BD do{
      xnew(ij) = (x(ij+(0,1)) + x(ij+(0,-1)) 
                   + x(ij+(1,0)) + x(ij+(-1,0))) / 4.0;

      xcheck(ij) = here.id;
    }

    delta = max reduce abs(xnew[ID] - x[ID]);
    x[ID] = xnew[ID];

    cnt += 1;
    if (verbose) {
      writeln("Iter: ", cnt, " (delta=", delta, ")\n");
      writeln(x);
      writeln("Jacobi mapping:");
      writeln(xcheck);
    }
  } while (delta > epsilon);
  return cnt;
}

proc gauss_seidel(D: domain(2), x: [D] real, epsilon: real) {
  const ID = D.expand(-1, -1); //domain for interior points
  var cnt = 0;                 //iteration counter

  //new array to check that mapping is correct
  var xcheck: [D] real = 8.8;

  //mapping   
  const BD = ID dmapped Block(boundingBox=ID);

  //boolean variables to check if delta < epsilon
  //and if any one delta is > epsilon.  If any is >
  //epsilon, need to keep calculating
  var flag1: bool = false;
  var flag2: bool = false;


  do {
    //initialize on each iteration
    flag1 = false;
    flag2 = false;

    forall ij in BD with (ref flag1, ref flag2) do{ 
      var temp: real = x(ij);

      x(ij) = (x(ij+(0,1)) + x(ij+(0,-1)) 
                 + x(ij+(1,0)) + x(ij+(-1,0))) / 4.0;

      xcheck(ij) = here.id;

      var my_delta = abs(x(ij) - temp);
      if (my_delta < epsilon){
        flag1 = true;
      }
      else{
        flag2 = true;
      }
    }

    cnt += 1;

    if (verbose) {
      writeln("Iter: ", cnt, "\n");
      writeln(x);
      writeln("Gauss-Seidel mapping:");
      writeln(xcheck);
    }

  } while (!(flag1 && !flag2));

  return cnt;
}

proc red_black(D: domain(2), x: [D] real, epsilon: real) {

  //domain for interior points
  const ID = D.expand(-1,-1); 

  //even domains
  const E1 = {2..n-2 by 2, 2..n-2 by 2};
  const E2 = {1..n-2 by 2, 1..n-2 by 2};

  //odd domains
  const O1 = {2..n-2 by 2, 1..n-2 by 2};
  const O2 = {1..n-2 by 2, 2..n-2 by 2};

  //counter
  var cnt = 0;

  //new array to check on mapping
  var xcheck: [D] real = 8.8;   

  //mapping evens and odds
  const BE1 = E1 dmapped Block(boundingBox=ID); //mapping
  const BE2 = E2 dmapped Block(boundingBox=ID); //mapping
  const BO1 = O1 dmapped Block(boundingBox=ID); //mapping
  const BO2 = O2 dmapped Block(boundingBox=ID); //mapping

  //flags to check if delta < epsilon and if any one delta > epsilon.
  //If any is > epsilon, need to keep calculating
  var flag1: bool;
  var flag2: bool;

  do {
    //initialize on each iteration
    flag1 = false;
    flag2 = false;

    forall ij in BE1 with (ref flag1, ref flag2) do {
      var temp: real = x(ij);

      x(ij) = (x(ij+(0,1)) + x(ij+(0,-1)) 
                 + x(ij+(1,0)) + x(ij+(-1,0))) / 4.0;

      xcheck(ij) = here.id;

      var my_delta = abs(x(ij) - temp);

      if (my_delta < epsilon){
        flag1 = true;
      }
      else{
        flag2 = true;
      }
    }

    forall ij in BE2 with (ref flag1, ref flag2) do {
      var temp: real = x(ij);

      x(ij) = (x(ij+(0,1)) + x(ij+(0,-1)) 
                 + x(ij+(1,0)) + x(ij+(-1,0))) / 4.0;

      xcheck(ij) = here.id;

      var my_delta = abs(x(ij) - temp);

      if (my_delta < epsilon){
        flag1 = true;
      }
      else{
        flag2 = true;
      }
    }

    forall ij in BO1 with (ref flag1, ref flag2) do {
      var temp: real = x(ij);

      x(ij) = (x(ij+(0,1)) + x(ij+(0,-1)) 
                 + x(ij+(1,0)) + x(ij+(-1,0))) / 4.0;

      xcheck(ij) = here.id;

      var my_delta = abs(x(ij) - temp);

      if (my_delta < epsilon){
        flag1 = true;
      }
      else{
        flag2 = true;
      }
    }

    forall ij in BO2 with (ref flag1, ref flag2) do {
      var temp: real = x(ij);

      x(ij) = (x(ij+(0,1)) + x(ij+(0,-1)) 
                 + x(ij+(1,0)) + x(ij+(-1,0))) / 4.0;

      xcheck(ij) = here.id;

      var my_delta = abs(x(ij) - temp);

      if (my_delta < epsilon){
        flag1 = true;
      }
      else{
        flag2 = true;
      }
    }
    
    cnt += 1;

    if (verbose) {
      writeln("Iter: ", cnt, "\n");
      writeln(x);
      writeln("Red-Black mapping:");
      writeln(xcheck);
    }

  } while (!(flag1 && !flag2)); 
  return cnt;
}



// Main routine.
//
proc main() {
  const D = {0..n-1, 0..n-1};   // domain including boundary points
  var a: [D] real = 0.0;	// mesh array
  a[n-1, 0..n-1] = 1.0;         // - setting boundary values
  a[0..n-1, n-1] = 1.0;
  var cnt = jacobi(D, a, epsilon);
  writeln("Jacobi:");
  writeln("Mesh size: ", n, " x ", n, ", epsilon=", epsilon, 
            ", total Jacobi iterations: ", cnt);
//  writeln(a);
  writeln("");

  //reset array
  a = 0.0;
  a[n-1, 0..n-1] = 1.0;
  a[0..n-1, n-1] = 1.0;
 
  //gauss-seidel method
  var gs_cnt = gauss_seidel(D, a, epsilon);
  writeln("Gauss-Seidel:");
  writeln("Mesh size: ", n, " x ", n, ", epslon: ", epsilon, 
          ", Gauss-Seidel iterations: ", gs_cnt);
//  writeln(a);
  writeln("");

  //reset array
  a = 0.0;
  a[n-1, 0..n-1] = 1.0;
  a[0..n-1, n-1] = 1.0;

  //red-black method
  var rb_cnt = red_black(D, a, epsilon);
  writeln("Red-Black:");
  writeln("Mesh size: ", n, " x ", n, ", epsilon: ", epsilon, 
          ", iterations: ", rb_cnt);
//  writeln(a);
  writeln("");

}
