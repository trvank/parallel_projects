//------------------------------------------------------------------------- 
// This is supporting software for CS415/515 Parallel Programming.
// Copyright (c) Portland State University
//------------------------------------------------------------------------- 

// Domain Map Examples
//
// Usage: ./domMap -nl <#locales>
// 
//

use BlockDist, CyclicDist;

config const n = 8;
const D = {1..n, 1..n};
const D2 = D.expand(-1, -1);

// A 2D Block-distributed domain BlockD and a Block-distributed 
// array BA declared over the domain.
//
const BlockD = D2 dmapped Block(boundingBox=D2);
var BA: [D] int = 0;

// A 2D Cyclic-distributed domain CyclicD and a Cyclic-distributed 
// array CA declared over the domain.
//
const CyclicD = D2 dmapped Cyclic(startIdx=D.low);
var CA: [D] int = 0;

// To illustrate how the index set is distributed across locales,
// we'll use forall loop to initialize each array element to the
// locale ID that stores that element.
//


forall e in BlockD do {
  BA(e) = here.id;
}

forall e in CyclicD do
  CA(e) = here.id;


// Output the arrays to visually see how the elements are
// partitioned across the locales.
//
writeln("Block Array:");
writeln(BA);

writeln("Cyclic Array:");
writeln(CA);
