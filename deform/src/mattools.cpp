#include "mattools.h"

int MatTools::cholmod_dense_mult(cholmod_dense* K, 
				 cholmod_dense* b, 
				 cholmod_dense* v) {
    
    char TRANSA = 'N'; // Implies that A and B are not transposed
    char TRANSB = 'N'; // Note that might be "T" because of alternative row/column ordering
    int M = K->nrow; // Rows of A
    int N = b->ncol; // Cols of B
    int _K = K->ncol; // Cols of A and Rows of B
    double alpha = 1.0; // A is scaled by this much
    double* A = (double*) K->x; // The matrix A
    int LDA = M; // The leading dimension of matrix A
    double* B = (double*) b->x; // The matrix B
    int LDB = _K; // Leading dimension of matrix B
    double beta = 0.0; // C is scaled by this much
    double* C = (double*) v->x; // The (output) matrix C
    int LDC = M; // Leading dimension of matrix C	
    //fprintf(stderr,"\ndgemm_(%c,%c,%d,%d,%d,%f,A,%d,B,%d,%f,C,%d)",
    //        TRANSA,TRANSB,M,N,_K,alpha,LDA,LDB,beta,LDC);
    dgemm_(&TRANSA, &TRANSB, &M, &N, &_K, &alpha, A, &LDA, B, &LDB, &beta, C, &LDC);    
    return 0;
}

int MatTools::cholmod_dense_solve(cholmod_dense* K, cholmod_dense* b) {
    // Note that the matrices passed will be mangled by this function!
    int N = K->nrow; // Order of A (number of linear equations in A)
    int NRHS = b->ncol; // Number of RHS's (i.e. number of solution columns)
    double *A = (double*) K->x; // Input matrix (to be factorised)
    int LDA = K->nrow; // Leading dimension of A
    int* IPIV = new int[N]; // Row swap indices (not used)
    double *B = (double*) b->x; // Answer vector / matrix
    int LDB = b->nrow; // Leading dimension of B
    int INFO = 0; // Output of this function
    dgesv_(&N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO);
    delete [] IPIV;
    return INFO;	   
}

/*
*/
int MatTools::cholmod_dense_solvels(cholmod_dense* K, cholmod_dense* b) {
    char TRANS = 'N'; // Whether the matrix is transposed
    int M = K->nrow; // Rows in matrix A
    int N = K->ncol; // Cols in matrix A
    int NRHS = b->ncol; // Number of columns of b
    double *A = (double*) K->x; // The matrix A
    int LDA = M; // Leading dimension of A
    double *B = (double*) b->x; // The matrix B
    int LDB = b->nrow; // The leading dimension of B
    double OPTWORK = 0; // Returned by a workspace query call
    int LWORK = -1; // Set to -1 for a workspace query, else use OPTWORK
    int INFO = 0; // Output of function (0 means success)
    // Perform a workspace size query
    dgels_(&TRANS, &M, &N, &NRHS, A, &LDA, B, &LDB, &OPTWORK, &LWORK, &INFO);
    if (INFO == 0) {
	double* WORK = new double[(int) OPTWORK]; // Allocate the optimal work size
	LWORK = (int) OPTWORK; // Set LWORK
	// Perform the actual LS solution
	dgels_(&TRANS, &M, &N, &NRHS, A, &LDA, B, &LDB, WORK, &LWORK, &INFO);
	delete [] WORK; // Clear away used memory
    }
    return INFO;
}


/** Invert the square matrix A of size N. Note that A will be modified with
 * this call. This function performs an LU decomposition with the LAPACK routine
 * DGETRF (modifying A into LU decomposition) and then performs the inverse
 * with DGETRI. Note that the contents of A cannot be trusted should this function
 * fail (failure is indicated by a return value != 0).
 * \param N The size of the square matrix A
 * \param A The square matrix to be inverted (MODIFIED)
 * \returns 0 on success, another value on failure (not sure how to interpret this value)
 */
int MatTools::invert(int N, double* A) {
    int LWORK = 4*N; // Define memory used by function
    int INFO = 0;       // Define the default info value
    int *IWORK = new int[N];
    double *WORK = new double[LWORK]; // Define a workspace for function
    int *IPIV = new int[N];

    // First compute the LU decomposition (stored in A) and the pivot vector
    // IPIV, which is needed for the inversion call
    // SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
    dgetrf_(&N, &N, A, &N, IPIV, &INFO);

    // Now determine if the condition number is acceptible
    double ANORM = norm(N, A);
    double RCOND = 0.0;
    dgecon_((char*) "1", &N, A, &N, &ANORM, &RCOND, WORK, IWORK, &INFO);

    if ((RCOND < EPSILON) || (INFO != 0) ) {
        // If things didn't go well, die gracefully
        // Clean up
        delete [] WORK;
        delete [] IPIV;
        delete [] IWORK;
        return -1;
    }

    // Now invert the matrix
    dgetri_(&N, A, &N, IPIV, WORK, &LWORK, &INFO );

    // A contains the inverted matrix, so we can simply go home
    //delete [] _A;
    delete [] WORK;
    delete [] IPIV;
    delete [] IWORK;
    return INFO;
}

/** Return the L1 norm of a matrix given its size
 */
double MatTools::norm(int N, double* A) {
    int i,j;
    double max = 0.0;
    double tot;
    for (i = 0; i < N; ++i) {
        tot = 0.0;
        for (j = 0; j < N; ++j) {
            tot += fabs(A[j*N+i]);
        }
        if (tot > max) max = tot;
    }
    return max;
}

/** Compute the determinant of a square input matrix of size N */
double MatTools::det(int N, double* A) {
    int INFO = 0;       // Define the default info value
    int *IPIV = new int[N];

    // First compute the LU decomposition (stored in A) and the pivot vector
    // IPIV, which is needed for the inversion call
    // SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
    dgetrf_(&N, &N, A, &N, IPIV, &INFO);

    // Correct the sign according to the number of row swaps (?)
    int i; double sign = 1.0;
    for (i = 0; i < N; ++i) if (IPIV[i]-1 != i) sign *= -1.0;

    delete [] IPIV;

    // The determinant is now just the product of the diagonal elements
    double _det = A[0];

    for (i = 1; i < N; ++i) _det *= A[i*N + i];
    return sign * _det;
}


