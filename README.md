# LinearAlgebra
Assignments for UE18MA251, Linear Algebra
Coded purely from scratch using Python 3.6.x. 

## Assignment - 1
Given a matrix A for the system of equations `Ax = b`
- Display, number of equations/rows, number of unknowns/columns, and rank of the matrix. (Use built in functions if required)
- Is the matrix singular (if it is a square matrix)?
- For a given `b`, will the system `Ax = b` have
  - Unique solution.
  - Infinitely many solutions
  - No Solution.

- If it has a unique solution display the solution.

## Assignment 2
- For a given matrix A which has a unique solution for Ax=b for any b,  
  test whether the LU factorization method is faster than the Gaussian Elimination method.

Run a test in this manner:
    -   First read a matrix A with full rank/full column rank.
    -   Generate a hundred random b vectors.
    -   Find a solution for each b.
    -   Note down how long this approach takes.

-   Factorize A to LU.
-   Generate a hundred random b vectors.
-   Find a solution for each b.
-   Record how long does this new approach take.

Check if the claim that the first approach takes `100*O(n3/3 )` and the second approach takes `O(n3/3 ) + 99*O(n2)` is true.
