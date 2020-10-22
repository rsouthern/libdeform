#ifndef _MATTOOLS_H
#define _MATTOOLS_H


#include <math.h>
#include <suitesparse/cholmod.h>

#ifndef EPSILON
#define EPSILON 0.000001
#endif

// Give a prototype for the fortran function dgemm_ for matrix
// multiplication (from the BLAS)
extern "C" {
    // Matrix multiplication
    void dgemm_(char*, char*, int*, int*, int*, double*,
		double*, int*, double*, int*, double*, 
		double*, int*);
    // Matrix solving
    void dgesv_(int*, int*, double*,
		int*, int*, double*,
		int*, int*);
    // Solve in the least squares sense
    void dgels_(char*, int*, int*, 
		int*, double*, int*, 
		double*, int*, double*, 
		int*, int*);

    // For the matrix inversion
    void dgetrf_(int*, int*, double*, int*, int*, int*);
    void dgetri_(int*, double*, int*, int*, double*, int*, int*);
    //SUBROUTINE DGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK,INFO )
    void dgecon_(char*, int*, double*, int*, double*, double*, double*, int*, int*);
};

namespace MatTools {
    // Useful cholmod routines
    int cholmod_dense_mult(cholmod_dense*, 
			   cholmod_dense*, 
			   cholmod_dense*);

    int cholmod_dense_solve(cholmod_dense*,
			    cholmod_dense*);

    int cholmod_dense_solvels(cholmod_dense*,
			      cholmod_dense*);

    /// Invert matrix A (size N)
    int invert(int N, double* A);

    /// Compute the L1 norm (max column sum) of matrix A (size N)
    double norm(int N, double* A);

    /// Compute the determinant of a square matrix by LU decomposition
    double det(int N, double* A);

};


#endif //_MATTOOLS_H


