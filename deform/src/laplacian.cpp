#include "laplacian.h"

extern "C" {
#include <svdlib.h>
}

/**
  * A constructor. Doesn't do anything yet.
  */
Laplacian::Laplacian() : dim(3) {
    init = false;

    // Retrieve the Choldmod workspace from the singleton class
    c = CholmodSharedSingleton::Instance()->getCholmodCommon();

    // Start cholmod with the current workspace
    X = NULL;
}

/**
  * A destructor. Should be virtual.
  */
Laplacian::~Laplacian() {
    if (init) clear();
    CholmodSharedSingleton::stop();

}

/**
 * Clear away contents of this struct so it can be used again. Note that
 * the constraints are NOT cleared unless explicitely requested.
 */
void Laplacian::clear(bool clearConstraints) {
    //qDebug("Laplacian::clear()");
    if (init) {
        cholmod_free_factor(&fac, c);
        cholmod_free_dense(&L_trans_B, c);
        cholmod_free_dense(&X, c);
        cholmod_free_dense(&B, c);
        cholmod_free_sparse(&L_trans, c);
        if (clearConstraints) {
            anchorList.clear();
            handleList.clear();
            handleMap.clear();
        }
        init = false;
    }

}

/**
  * A couple of fairly trivial functions for setting and
  * unsetting constraints.
  */
void Laplacian::addHandle(int i) {
    if (init) return;
    list<int>::iterator it = find(handleList.begin(), handleList.end(), i);
    if (it == handleList.end()) handleList.push_back(i);
}

void Laplacian::addAnchor(int i) {
    if (init) return;
    list<int>::iterator it = find(anchorList.begin(), anchorList.end(), i);
    if (it == anchorList.end()) anchorList.push_back(i);
}

/**
 * Retrieve a pointer to the deformed coordinates.
 */
cholmod_dense *Laplacian::getX() {
    if (!init) {
        printf("Laplacian::getX() - trying to retrieve an uninitialised value of X!\n");
        exit(0);
    }
    return X;
}

/**
  * Set up the relevant laplacian matrices to be ready for deformation.
  * The heavy, labour intensive stuff should be done here.
  */
void Laplacian::initialise(Mesh *m) {
    //clock_t start, cuda_time, cpu_time;
    //cuda_time = cpu_time = 0;
    //start = clock();
    //cuda_time += clock() - start;

    //start = clock();
    //cpu_time += clock() - start;




    if ((init)||(anchorList.empty() && handleList.empty())) return;


    // Set an initial version of X
    X = cholmod_copy_dense(m->getX(), c);



    float SECS_PER_CLOCK = 1.0 / (float) CLOCKS_PER_SEC;


    // Construct the LHS matrix L and add constraints
    clock_t start = clock();
    cholmod_sparse *L = constructLHS(m);
    fprintf(stderr, "\ - nconstructLHS() - ELAPSED TIME %fs", (clock()-start)*SECS_PER_CLOCK);


    // Construct the RHS (B for simple Laplacian)
    start = clock();
    constructRHS(L);
    fprintf(stderr, "\n - constructRHS() - ELAPSED TIME %fs", (clock()-start)*SECS_PER_CLOCK);

    // Factorise and ready the problem for solving
    start = clock();
    prepareProblem(L);
    fprintf(stderr, "\n - prepareProblem() - ELAPSED TIME %fs", (clock() - start)*SECS_PER_CLOCK);

    // Set this to be initialised (so we can deform it)
    init = true;

    // Clear away allocated memory
    cholmod_free_sparse(&L, c);
}

/**
  * Construct the LHS matrix. In the simple Laplacian case this is laplacian matrix L.
  */
cholmod_sparse *Laplacian::constructLHS(Mesh *m) {
    // Store some valuable variables
    vertNum = m->n_vertices();

    // Clear away our anchor rows
    handleMap.clear();

    // This is pretty inefficient, as we don't know exactly how big the
    // one rings are going to be around each vertex so we set nzmax to
    // the size of the matrix for the time being. We should count the
    // number of entries and then doublelocate the triplet at the end.
    cholmod_triplet *L_trip =
            cholmod_allocate_triplet(
                    vertNum+anchorList.size() + handleList.size(), //size_t nrow, /* # of rows of T */
                    vertNum, //size_t ncol, /* # of columns of T */
                    vertNum*dim*dim+anchorList.size()+handleList.size(), // size_t nzmax, /* max # of nonzeros of T */
                    0, //int stype, /* stype of T (0 is unsymmetric) */
                    CHOLMOD_REAL, //int xtype,/* CHOLMOD_PATTERN, _REAL, _COMPLEX, or _ZOMPLEX */
                    c);

    if (L_trip != NULL) {
        //cholmod_print_triplet(L_trip, "L_trip", &c);
    } else {
        printf("\nFailure to allocate L_trip(%d,%d), status=%d!", (int) vertNum+((int) handleList.size() + anchorList.size()), vertNum, c->status);
        exit(0);

    }

    // Calculate our differential coordinate matrices L and B
    diffMatrix(m, L_trip);

    int* L_tripi = (int*) L_trip->i;
    int* L_tripj = (int*) L_trip->j;
    double* L_tripx = (double*) L_trip->x;

    // Adding our constraints to the end of our matrix L
    list<int>::iterator it;
    int cnt;

    // Add our anchors first
    for (it = anchorList.begin(), cnt = vertNum; it != anchorList.end(); ++it, ++cnt) {
        // Do a little dynamic resizing if necessary
        if (L_trip->nnz+1 >= L_trip->nzmax) {
            cholmod_reallocate_triplet(L_trip->nnz + 64, L_trip, c);
        }
        // Append our anchor
        L_tripi[L_trip->nnz] = cnt;
        L_tripj[L_trip->nnz] = *it;
        L_tripx[L_trip->nnz] = 1.0;
        L_trip->nnz++;
    }

    // Now add our handles
    for (it = handleList.begin(); it != handleList.end(); ++it, ++cnt) {
        // Do a little dynamic resizing if necessary
        if (L_trip->nnz+1 >= L_trip->nzmax) {
            cholmod_reallocate_triplet(L_trip->nnz + 64, L_trip, c);
        }
        // Append our anchor
        L_tripi[L_trip->nnz] = cnt;
        L_tripj[L_trip->nnz] = *it;
        L_tripx[L_trip->nnz] = 1.0;
        L_trip->nnz++;

        // Set our map value
        handleMap[*it] = cnt;
    }

    // Make L_trip into a sparse
    cholmod_sparse *L =
            cholmod_triplet_to_sparse(
                    L_trip,               /* matrix to copy */
                    L_trip->nnz,          /* allocate at least this much space in output matrix */
                    c);
    cholmod_free_triplet(&L_trip, c);
    return L;
}

/**
  * Check the left hand side of the problem by performing sparse SVD and checking that the singular values are all
  * less than one.
  */
bool Laplacian::isMatValid(cholmod_sparse *L) {
    // Copy the cholmod_sparse structure into an svdlib structure
    smat Ls;
    int i;

    // These are easy to do
    Ls.rows = L->nrow;
    Ls.cols = L->ncol;
    Ls.vals = L->nzmax;

    // Need to add one to the column pointer (I think)
    Ls.pointr = (long*) malloc(sizeof(long)*L->ncol+1);
    for (i=0; i<=L->ncol; ++i) Ls.pointr[i] = (long) (((int*) L->p)[i]);

    // Need to typecast the data in the L->i field
    Ls.rowind = (long*) malloc(sizeof(long)*L->nzmax);
    for (i=0; i<L->nzmax; ++i) Ls.rowind[i] = (long) (((int*)L->i)[i]);

    // Just do a raw copy of the x data
    Ls.value = (double*) malloc(sizeof(double)*L->nzmax);
    memcpy(Ls.value, L->x, sizeof(double)*L->nzmax);

    // Set up the parameters for SVD
    double las2end[2] = {-1.0e-30, 1.0e-30}; // Values that are close enough to zero to be zero
    double kappa = 1e-6; // Tolerance for algorithm
    int dimensions = (Ls.cols>Ls.rows)?Ls.rows:Ls.cols; // Don't know
    int iterations = 5; // Number of iterations

    SVDRec R = NULL;
    R = svdLAS2(&Ls, dimensions, iterations, las2end, kappa);

    if (R != NULL) {
        // Solution was successful: now we analyse the singular values
        fprintf(stderr, "\nSingular Values: ");
        for (i=0; i<R->d; ++i) {
            fprintf(stderr,"%f, ", R->S[i]);
        }
        svdFreeSVDRec(R);
    }

    // Clear away memory
    free(Ls.pointr);
    free(Ls.rowind);
    free(Ls.value);

    return true;
}

/**
  * Compute B = L * X_0
  */
void Laplacian::constructRHS(cholmod_sparse *L) {
    // Just allocate a big chunk for our dense matrix B
    B = cholmod_allocate_dense(
            vertNum+handleList.size()+anchorList.size(), // size_t nrow, /* # of rows of matrix */
            dim,    // size_t ncol, /* # of columns of matrix */
            vertNum+handleList.size()+anchorList.size(), // size_t d, /* leading dimension */
            CHOLMOD_REAL,  // int xtype, /* CHOLMOD_REAL, _COMPLEX, or _ZOMPLEX */
            c);
    double alpha[2] = {1.0, 0.0};
    double beta[2] = {0.0, 0.0};
    cholmod_sdmult(L, 0, alpha, beta, X, B, c);
}

/**
  * Prepare the problem using the preconstructed LHS matrix L. Note that B is stored in this
  * structure, but L is not.
  */
void Laplacian::prepareProblem(cholmod_sparse *L) {
    // Make a transposed version of L
    L_trans = cholmod_allocate_sparse(
            L->ncol, /* # of rows of A */
            L->nrow, /* # of columns of A */
            L->nzmax,  /* max # of nonzeros of A */
            FALSE,  /* TRUE if columns of A sorted, FALSE otherwise */
            FALSE,  /* TRUE if A will be packed, FALSE otherwise */
            0,       /* stype of A */
            CHOLMOD_REAL, /* CHOLMOD_PATTERN, _REAL, _COMPLEX, or _ZOMPLEX */
            c);

    cholmod_transpose_unsym(
            L, /* matrix to transpose */
            1, /* 0: pattern, 1: array transpose, 2: conj. transpose */
            NULL, /* size nrow, if present (can be NULL) */
            NULL, /* subset of 0:(A->ncol)-1 */
            0,    /* size of fset */
            L_trans, /* F = A’, A(:,f)’, or A(p,f)’ */
            c);

    // Create the matrix L'L
    cholmod_sparse *L_trans_L = cholmod_ssmult(
            L_trans, /* left matrix to multiply */
            L,       /* right matrix to multiply */
            0,       /* requested stype of C */
            TRUE,   /* TRUE: do numerical values, FALSE: pattern only */
            TRUE,   /* if TRUE then return C with sorted columns */
            c);

    // Verify the LHS of the system by performing SVD on the sparse matrix
    // if (!isMatValid(L_trans_L)) {}

    // Analyse, factorise  and solve L'L
    L_trans_L->stype = 1; // Convert L_trans_L to be upper symmetric
    fac = cholmod_analyze(L_trans_L, /* matrix to order and analyze */
                          c);
    cholmod_factorize(
            L_trans_L, /* matrix to factorize */
            fac,       /* resulting factorization */
            c);            

    // Create L'B
    L_trans_B = cholmod_allocate_dense(
            L_trans->nrow, // size_t nrow, /* # of rows of matrix */
            B->ncol,    // size_t ncol, /* # of columns of matrix */
            L_trans->nrow, // size_t d, /* leading dimension */
            CHOLMOD_REAL,  // int xtype, /* CHOLMOD_REAL, _COMPLEX, or _ZOMPLEX */
            c);

    // Clear up memory
    cholmod_free_sparse(&L_trans_L, c);
}

/**
 * Construct a matrix of differential coordinates for this complex. Note that
 * it doesn't have to be a manifold - it could be a full space complex.
 * The form of our system is Ax=B where A is a n*n matrix and B is a n*m matrix,
 * where n is the number of vertices and m is the dimension of the vertices.
 */
void Laplacian::diffMatrix(Mesh *m, cholmod_triplet *A) {
    if (init) return;

    // Finding the differential coordinates for each vertex
    int i,idx;

    int* Ai = (int*) A->i;
    int* Aj = (int*) A->j;
    double* Ax = (double*) A->x;

    uint oRsize;

    Mesh::VertexIter v_it;
    Mesh::VertexVertexIter vv_it;

    for (v_it = m->vertices_begin(), i=0; v_it != m->vertices_end(); ++v_it, ++i) {
        // Count the size of the one ring (is there a better way to do this?
        oRsize = 0;
        for (vv_it = m->vv_iter(v_it.handle()); vv_it; ++vv_it, ++oRsize);
        double _w = 1.0 / ((double) oRsize);

        // Make sure there is space in our triplet
        if ((A->nnz + oRsize + 1) >= A->nzmax) {
            cholmod_reallocate_triplet(A->nnz+oRsize+2, A, c);
        }

        // Iterate over all one ring vertices
        for (vv_it = m->vv_iter(v_it.handle()), idx=0; vv_it; ++vv_it,++idx) {
            int v_idx = vv_it.handle().idx();
            // Storing the fraction into the matrix A
            Ai[A->nnz] = i;
            Aj[A->nnz] = v_idx;
            Ax[A->nnz] = _w;
            A->nnz++;
        }

        // Add our i entry in this row
        Ai[A->nnz] = i;
        Aj[A->nnz] = i;
        Ax[A->nnz] = -1.0;
        A->nnz++;
    }
}

/**
  *
  */
void Laplacian::transformHandles(cholmod_dense *M) {
    if (M != NULL) {   
        // Allocate a chunk for our anchors
        cholmod_dense *_B = cholmod_allocate_dense(
                M->nrow,                // size_t nrow,
                handleList.size(),      // size_t ncol,
                M->nrow,                // size_t d,
                CHOLMOD_REAL,           // int xtype, /* CHOLMOD_REAL, _COMPLEX, or _ZOMPLEX */
                c);

        // Allocate a chunk for the deformed anchors
        cholmod_dense *M_B = cholmod_allocate_dense(
                M->nrow,                // size_t nrow, /* # of rows of matrix */
                handleList.size(),      // size_t ncol, /* # of columns of matrix */
                M->nrow,                // size_t d, /* leading dimension */
                CHOLMOD_REAL,           // int xtype, /* CHOLMOD_REAL, _COMPLEX, or _ZOMPLEX */
                c);


        double* _Bx = (double*) _B->x;
        double* Bx = (double*) B->x;

        list<int>::iterator it;
        int j, cnt;

        // Create the matrix of anchors
        for (it = handleList.begin(),cnt=0; it != handleList.end(); ++it,++cnt) {
            for (j = 0; j < dim; ++j) {
                _Bx[j + cnt*_B->nrow] = Bx[handleMap[*it]+j*B->nrow];
            }
            _Bx[dim + cnt*_B->nrow] = 1.0;
        }

        MatTools::cholmod_dense_mult(M, _B, M_B);
        double* M_Bx = (double*) M_B->x; // The (output) matrix C

        // Copy the result back into B
        for (it = handleList.begin(), cnt=0; it != handleList.end(); ++it,++cnt) {
            for (j = 0; j < dim; ++j) {
                if (M_Bx[dim + cnt*M_B->nrow] == 0.0) {
                    fprintf(stderr, "\nLaplacian::transformAnchors() - transformed points have zero denominator!");
                } else {                    
                    Bx[handleMap[*it]+j*B->nrow] = M_Bx[j + cnt*M_B->nrow] / M_Bx[dim + cnt*M_B->nrow];
                }
            }
        }

        //cholmod_print_dense(M, "M", c);
        //cholmod_print_dense(_B, "_B", c);
        //cholmod_print_dense(M_B, "M_B", c);

        // Clear away allocated memory
        cholmod_free_dense(&_B, c);
        cholmod_free_dense(&M_B, c);
    }
}

/**
  * Recompute the laplacian solution L'LX=L'B. This implies that B has already been computed
  */
void Laplacian::update() {
    if (!init) return;

    // Multiply L' (sparse) by B
    double alpha[2] = {1.0, 0.0};
    double beta[2] = {0.0, 0.0};
    cholmod_sdmult(
            L_trans,     /* sparse matrix to multiply */
            0,           /* use A if 0, or A’ if 1, or A.’ if -1 */
            alpha,       /* scale factor for A */
            beta,        /* scale factor for Y */
            B,           /* dense matrix to multiply */
            L_trans_B,   /* resulting dense matrix */
            c);

    if (X != NULL) cholmod_free_dense(&X, c);
    X = cholmod_solve(
            0,         /* system to solve */
            fac,       /* factorization to use */
            L_trans_B, /* right-hand-side */
            c);

    //cholmod_print_dense(X, (char*) "X", c);
}



/**
  * Print this structure to the screen or somewhere else appropriate
  */
void Laplacian::print() {
    FILE *fid = stderr;
    if (init) {
        cholmod_print_dense(X, "X", c);
        cholmod_print_dense(B, "B", c);
        cholmod_print_sparse(L_trans, "L_trans", c);
        cholmod_print_dense(L_trans_B, "L_trans_B", c);
        cholmod_print_factor(fac, "fac", c);
    }
}


void Laplacian::saveConstraints(FILE *fid) {
    fprintf(fid, "%d %d\n", (int) anchorList.size(), (int) handleList.size());
    list<int>::iterator it;
    for (it = anchorList.begin(); it != anchorList.end(); ++it) {
        fprintf(fid, "%d\n", *it);
    }
    for (it = handleList.begin(); it != handleList.end(); ++it) {
        fprintf(fid, "%d\n", *it);
    }
}

void Laplacian::loadConstraints(FILE *fid) {
    if (!init) {
        int nA, nH, v, i;

        fscanf(fid, "%d %d", &nA, &nH);
        for (i=0; i<nA; ++i) {
            fscanf(fid, "%d", &v);
            anchorList.push_back(v);
        }
        for (i=0; i<nH; ++i) {
            fscanf(fid, "%d", &v);
            handleList.push_back(v);
        }
    }
    //cerr << "\nLaplacian::loadConstraints() - anchors=["<<anchorList<<"], handles=["<<handleList<<"]";
}
