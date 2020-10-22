#include "laplacian_DUAL.h"
#include <omp.h>
#include "cusp_tools.cuh"

/**
 *
 */
Laplacian_DUAL::Laplacian_DUAL() : Laplacian() {
}

Laplacian_DUAL::~Laplacian_DUAL() {
}


/**
 * Perform a clearup. Call the parent class too while we're there.
 */
void Laplacian_DUAL::clear(bool _c) {
    if (init) {
        cholmod_free_dense(&DX,c);
        Laplacian::clear(_c);
    }
}

/**
 * Set up the relevant laplacian matrices to be ready for deformation.
 * The heavy, labour intensive stuff should be done here.
 */
void Laplacian_DUAL::initialise(Mesh *m) {
    if ((init)||(anchorList.empty() && handleList.empty())) return;

    // Construct our dual mesh using the primal mesh
    dm.construct(m);

    // Make a local copy of the dual vertices
    DX = cholmod_copy_dense(dm.getX(), c);

    // Call the parent class initialise function
    Laplacian::initialise(m);
}

cholmod_sparse *Laplacian_DUAL::constructLHS(Mesh *m) {
    // Clear away our anchor rows
    handleMap.clear();

    // We construct the dual Laplacian first, consisting of L = [\hat{L}*D; I_c]
    cholmod_triplet *L_trip = cholmod_allocate_triplet(
                    dm.n_vertices()+anchorList.size()+handleList.size(), //size_t nrow, /* # of rows of T */
                    dm.getD()->ncol,           //size_t ncol, /* # of columns of T */
                    (dm.n_vertices()*20)+anchorList.size()+handleList.size(), // size_t nzmax, /* max # of nonzeros of T */
                    0, //int stype, /* stype of T (0 is unsymmetric) */
                    CHOLMOD_REAL, //int xtype,/* CHOLMOD_PATTERN, _REAL, _COMPLEX, or _ZOMPLEX */
                    c);

    // Calculate our differential coordinate matrices L and B
    diffMatrix(m, L_trip);

    int* L_tripi = (int*) L_trip->i;
    int* L_tripj = (int*) L_trip->j;
    double* L_tripx = (double*) L_trip->x;

    // Adding our constraints to the end of our matrix L
    list<int>::iterator it;
    int cnt;

    // Add our anchors first
    for (it = anchorList.begin(), cnt = dm.n_vertices(); it != anchorList.end(); ++it, ++cnt) {
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
                    L_trip, /* matrix to copy */
                    L_trip->nnz,          /* allocate at least this much space in output matrix */
                    c);
    cholmod_free_triplet(&L_trip, c);
    return L;
}

/**
  * Compute B = L * X_0
  */
void Laplacian_DUAL::constructRHS(cholmod_sparse *L) {
    // Compute the answer vector B
    // Just allocate a big chunk for our dense matrix B
    B = cholmod_allocate_dense(
            dm.n_vertices()+anchorList.size()+handleList.size(), // size_t nrow, /* # of rows of matrix */
            dim,    // size_t ncol, /* # of columns of matrix */
            dm.n_vertices()+anchorList.size()+handleList.size(), // size_t d, /* leading dimension */
            CHOLMOD_REAL,  // int xtype, /* CHOLMOD_REAL, _COMPLEX, or _ZOMPLEX */
            c);
    double alpha[2] = {1.0, 0.0};
    double beta[2] = {0.0, 0.0};
    cholmod_sdmult(L, 0, alpha, beta, X, B, c);
}

/**
 * Construct a matrix of differential coordinates for this complex. Note that
 * it doesn't have to be a manifold - it could be a full space complex.
 * The form of our system is Ax=B where A is a n*n matrix and B is a n*m matrix,
 * where n is the number of vertices and m is the dimension of the vertices.
 */
void Laplacian_DUAL::diffMatrix(Mesh *m, cholmod_triplet *A) {
    if (init) return;

    // Finding the differential coordinates for each vertex
    int i,idx;

    // Store the temporary dual laplacian matrix here
    cholmod_triplet *_L_trip = cholmod_allocate_triplet(
            dm.n_faces(), //size_t nrow,
            dm.n_faces(), //size_t ncol,
            dm.n_faces()*(dim+1), // size_t nzmax,
            0, //int stype,
            CHOLMOD_REAL, //int xtype,
            c);

    // Our random access into the dual matrix
    int* _Li = (int*) _L_trip->i;
    int* _Lj = (int*) _L_trip->j;
    double* _Lx = (double*) _L_trip->x;

    // Will store the vector of the iterative corrections d_i
    correction.resize(dm.n_faces());
    faceArea.resize(dm.n_faces());

    // We must now find the correction matrix and the laplacian
    vector<double> w_ij(dim);

    // Loop over all our dual facets
    Mesh::FaceIter f_it;
    Mesh::VertexIter v_it;
    Mesh::FaceVertexIter fv_it;

    for (f_it = dm.faces_begin(), v_it = dm.vertices_begin(), idx=0; f_it != dm.faces_end(); ++f_it, ++idx, ++v_it) {
        // Find the weights of the given simplex and its parent vertex
        findWeights(v_it.handle(), f_it.handle(), &w_ij, &(correction[idx]), &(faceArea[idx]));

        // Assign the weights to our dual laplacian        
        for (fv_it = dm.fv_iter(f_it.handle()),i=0;fv_it;++fv_it,++i) {
            _Li[_L_trip->nnz] = idx;
            _Lj[_L_trip->nnz] = fv_it.handle().idx();
            _Lx[_L_trip->nnz] = w_ij[i];
            _L_trip->nnz++;
        }

        // Let us hope that our barycentric coordinates sum to 1!
        _Li[_L_trip->nnz] = idx;
        _Lj[_L_trip->nnz] = idx;
        _Lx[_L_trip->nnz] = -1.0;
        _L_trip->nnz++;
    }
    //cholmod_print_triplet(_L_trip, "L_trip", c);


    // Now compute our true matrix A from the _L * D
    cholmod_sparse *_L = cholmod_triplet_to_sparse(
            _L_trip, 		// matrix to copy
            _L_trip->nnz,    // allocate at least this much space in output matrix
            c);
    cholmod_free_triplet(&_L_trip, c);

    cholmod_sparse *A_sparse = cholmod_ssmult(
            _L, // left matrix to multiply
            dm.getD(),// right matrix to multiply
            0,       // requested stype of C
            TRUE,   // TRUE: do numerical values, FALSE: pattern only
            TRUE,   // if TRUE then return C with sorted columns
            c);
    cholmod_free_sparse(&_L, c);
    cholmod_triplet *_A = cholmod_sparse_to_triplet(A_sparse, c);
    cholmod_free_sparse(&A_sparse, c);

    int* _Ai = (int*) _A->i;
    int* _Aj = (int*) _A->j;
    double* _Ax = (double*) _A->x;

    int* Ai = (int*) A->i;
    int* Aj = (int*) A->j;
    double* Ax = (double*) A->x;

    // Do a deep copy of the contents of _A into A (is there a faster way to do this?)
    memcpy(Ai, _Ai, sizeof(int)*_A->nnz);
    memcpy(Aj, _Aj, sizeof(int)*_A->nnz);
    memcpy(Ax, _Ax, sizeof(double)*_A->nnz);
    A->nnz = _A->nnz;

    // Clean up memory
    cholmod_free_triplet(&_A, c);
}

/**
 * This function finds the weight of a given vertex in a given simplex by finding
 * the barycentric coordinates of the vertex projected into the plane of the
 * simplex. The function will return two things - a vector containing the
 * barycentric coordinates, and a vector containing the value of -d_i n_i which
 * is the corrective vector used to construct the differential vector.
 * \return h_i, the magnitude of the vertex displacement.
 */
void Laplacian_DUAL::findWeights(VHandle vh, FHandle fh, vector<double> *w_ij, double *c, double *a) {
    vector<double> normal(3,0);
    (*a) = findNormal(fh, &normal);
    int i;
    int N = dim;

    // If this t is valid, we need to compute the barycentric coordinates
    double *A = new double[N*N];
    double *b = new double[N];

    // Fill up A with the vectors of this simplex
    // The formulation of A comes from
    // [ (b-a) (c-a) -d ] * [beta; alpha; t] = (p-a)
    Mesh::FaceVertexIter fv_it = dm.fv_iter(fh);
    Mesh::Point v0 = dm.point(fv_it.handle());

    // Set the first two columns of A to (b-a) and (c-a) respectively
    //cerr << "\nV0=[" << v0 << "]';";
    ++fv_it;
    for (i=0; i < 3; ++i) {
        A[i+0*N] = dm.point(fv_it.handle())[i] - v0[i];
    }
    //cerr << "\nV1=[" << dm.point(fv_it.handle())<< "]';";
    ++fv_it;
    for (i=0; i < 3; ++i) {
        A[i+1*N] = dm.point(fv_it.handle())[i] - v0[i];
    }
    //cerr << "\nV2=[" << dm.point(fv_it.handle())<< "]';";

    //cerr << "\nPt=["<<dm.point(vh)<< "]';";
    // Set the last column of A to -d, the normal
    for (i=0; i < 3; ++i) {
        A[i+2*N] = -normal[i];
    }
    // Set the answer vector B to (p-a)
    for (i=0; i<3; ++i) {
        b[i] = dm.point(vh)[i] - v0[i];
    }

    //fprintf(stderr,"\nA=[%f,%f,%f;%f,%f,%f;%f,%f,%f]'",A[0],A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8],A[9]);
    //fprintf(stderr,"\nb=[%f;%f;%f]",b[0],b[1],b[2]);

    // Solve the system for our barycentric coordinates
    int  *IPIV = new int[N];
    int NRHS = 1;
    int INFO = 0;
    dgesv_(&N, &NRHS, A, &N, IPIV, b, &N, &INFO);

    //fprintf(stderr, "\nres=[%f;%f;%f]",b[0],b[1],b[2]);

    // If INFO is nonzero, then there was a failure - not sure why this could
    // happen. The vertex + normal must intersect the plane somewhere. If it
    // failed, we simply don't make a baryCoord record and return.
    double h_i = 0.0;

    if (INFO == 0) {
        h_i = - b[2]; // Is this correct?
        // Now we find the final barycentric coordinate by totalling the others
        w_ij->at(0) = 1.0 - (b[0] + b[1]);
        w_ij->at(1) = b[0];
        w_ij->at(2) = b[1];
    } else {
        fprintf(stderr,"\nLaplacian::findWeights() - barycentric coordinate computation failed!");
        exit(0);
    }

    (*c) = h_i;   

    delete [] IPIV;
    delete [] b;
    delete [] A;
}

/**
 * Here we apply the iterative post-process which will allow our mesh to converge to
 * the minimum energy mesh. The process is not exactly the same as the original
 * dual laplacian paper -- we just use an absolute distance between the vertices
 * to decide on convergence, rather than a ratio of the distances between steps.
 * This should probably be changed at some point.
 * \see needRestore(), updateDualVerts(), correct(), deform()
 *
 */
void Laplacian_DUAL::update() {
    // Sanity check
    if (!init) return;

    // Perform one loop, storing the current vertex positions    
    cholmod_dense *_X = cholmod_copy_dense(X, c);

    updateDX();             // Update the dual vertices
    correct();              // Update the detail vectors according to local rotations
    Laplacian::update();
    double eps = epsilon(_X, X);
    int it = 0;
    int max_it = 100;
    double min_eps = 0.00001;
    fprintf(stderr,"\n - %d: epsilon(_X,X) = %f",it, eps);
    //char fname[256];

    while ((eps > min_eps) && (it < max_it)) {
        cholmod_copy_dense2(X, _X, c);
        updateDX();         // Update the dual vertices
        correct();          // update the detail vectors according to local rotations

        //dm.setX(DX);
        //sprintf(fname, "data/deform_%d.obj",it);
        //dm.save(fname);

        Laplacian::update();    // Call the parent update function to solve for X
        eps = epsilon(_X, X);   // Compute the error
        ++it;       
        fprintf(stderr,"\n - %d: epsilon(_X,X) = %f",it, eps);
    }
}

void Laplacian_DUAL::update_cuda() {
    // The update function may be called just to create an initial X
    if (X == NULL) {
        Laplacian::update();
        return;
    }

    // This is going to crash!
    //test_cusp(fac);

    // First initialise cuda to use our data
    CudaDual::init_dual_cu(DX->nrow, B->nrow, dm.n_faces(), dm.getTri(), &(correction[0]), &(faceArea[0]));

    // Perform one loop, storing the current vertex positions
    cholmod_dense *_X = cholmod_copy_dense(X, c);
    updateDX();         // Update the dual vertices

    // Update the detail vectors according to local rotations
    CudaDual::correct_dual_cu((double*)DX->x, (double*)B->x);

    Laplacian::update();
    double eps = epsilon(_X, X);

    int it = 0;
    int max_it = 100;
    double min_eps = 0.0001;
    fprintf(stderr,"\n - %d: epsilon(_X,X) = %f",it, eps);
    //char fname[256];

    //clock_t start, cuda_time, cpu_time;
    //cuda_time = cpu_time = 0;

    while ((eps > min_eps) && (it < max_it)) {
        cholmod_copy_dense2(X, _X, c);
        updateDX();         // Update the dual vertices
        //start = clock();
        CudaDual::correct_dual_cu((double*)DX->x, (double*)B->x);
        //cuda_time += clock() - start;

        //start = clock();
        //correct();
        //cpu_time += clock() - start;

        //cholmod_print_dense(B,"B",c);

        //dm.setX(DX);
        //sprintf(fname, "debug_%d.obj",it);
        //dm.save(fname);

        Laplacian::update();    // Call the parent update function to solve for X
        eps = epsilon(_X, X);   // Compute the error
        ++it;
        fprintf(stderr,"\n - %d: epsilon(_X,X) = %f",it, eps);
    }

    CudaDual::shutdown_dual_cu();


    //float SECS_PER_CLOCK = 1.0 / (float) CLOCKS_PER_SEC;
    //fprintf(stderr,"\nCUDA %fs, CPU %fs", ((float)cuda_time)*SECS_PER_CLOCK, ((float)cpu_time)*SECS_PER_CLOCK);
}

void Laplacian_DUAL::updateDX() {
    // Update the dual vertices - simply multiply the dual matrix operator D our current vertices X
    double alpha[2] = {1.0, 0.0};
    double beta[2] = {0.0, 0.0};

    cholmod_sdmult(
            dm.getD(),   /* sparse matrix to multiply */
            0,           /* use A if 0, or A’ if 1, or A.’ if -1 */
            alpha,       /* scale factor for A */
            beta,        /* scale factor for Y */
            X,           /* dense matrix to multiply */
            DX,          /* resulting dense matrix */
            c);
}

/**
 * Correct the matrix B by assigning to it the value -h_i * n_i
 */
void Laplacian_DUAL::correct() {
    // Now we have the new vertices, we compute the normal for each face of
    // the dual complex, and update the matrix B with the result    
    double *Bx = (double*) B->x;
    int i,j;
    Mesh::FaceIter f_it;
    vector<double> normal(3,0);
    double a;
    double scale;

    for (f_it = dm.faces_begin(), i=0; f_it != dm.faces_end(); ++f_it, ++i) {
        a = findNormal(f_it.handle(), &normal);
        scale = sqrt(a/faceArea[i]);
        for (j = 0; j < (int) B->ncol; ++j) {
            Bx[i + j*B->nrow ] = - scale * correction[i] * normal[j];
        }
    }
}

/**
  * Our dual mesh does not hold the current vertex data - it is stored in X. Accordingly
  * we must update the face normal using the updated data.
  */
double Laplacian_DUAL::findNormal(Mesh::FaceHandle fh, vector<double> *normal) {
    // Find the index of the vertices in this face
    vector<double> p0(3,double(0)), v1(3,double(0)), v2(3,double(0)), res(3,double(0));
    vector<int> indices(3,double(0));
    Mesh::FaceVertexIter fv_it;
    int i,j,idx;
    for (fv_it = dm.fv_iter(fh),i=0; fv_it; ++fv_it,++i) {
        idx = fv_it.handle().idx();
        indices[i] = idx;
        if (i==0) {
            // set our first point
            for (j=0; j<3; ++j)
                p0[j] = ((double*)DX->x)[idx + j*DX->nrow];
        } else if (i==1) {
            // Create our vector
            for (j=0; j<3; ++j)
                v1[j] = ((double*) DX->x)[idx + j*DX->nrow] - p0[j];
        } else {
            for (j=0; j<3; ++j)
                v2[j] = ((double*) DX->x)[idx + j*DX->nrow] - p0[j];
        }
    }

    // Now compute the cross product
    CROSS3(res, v1, v2);

    // Normalise the result
    double tot=res[0]*res[0]+res[1]*res[1]+res[2]*res[2];
    double inv_tot = 1.0 / sqrt(tot);
    for (i=0; i<3; ++i) {
        (*normal)[i] = res[i] * inv_tot;
    }

    return sqrt(tot);
}

/**
 * Compute the maximum edge length as a norm of convergence. Assume A and B
 * are same dimensions!
 */
double Laplacian_DUAL::epsilon(cholmod_dense *first, cholmod_dense *second) {
    double max = 0.0;
    double dist;

    double *Ax = (double*)first->x;
    double *Bx = (double*)second->x;

    int i,j;
    for (i = 0; i < (int) first->nrow; ++i) {
        dist = 0.0;
        for (j = 0; j < (int) first->ncol; ++j) {
            dist += (Ax[j * first->nrow + i] - Bx[j * first->nrow + i]) * (Ax[j * first->nrow + i] - Bx[j * first->nrow + i]);
        }
        if ((i == 0) || (dist > max)) {
            max = dist;
        }
    }
    return max;
    //return sqrt(max);
}
