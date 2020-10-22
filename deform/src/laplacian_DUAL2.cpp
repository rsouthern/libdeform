#include "laplacian_DUAL2.h"

// Define this to avoid CGAL from bombing out when the system is not positive definite
// Also apparently this also improves performance
#define CGAL_QP_NO_ASSERTIONS

#include <CGAL/basic.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

// choose exact integral type
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpzf.h>
typedef CGAL::Gmpzf ET;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#endif

typedef CGAL::Const_oneset_iterator<CGAL::Comparison_result> Relation;

typedef CGAL::Quadratic_program_from_iterators
<ET**,                                                // for A
 ET*,                                                 // for b
 Relation,                                                // for r
 bool*,                                                   // for fl
 ET*,                                                 // for l
 bool*,                                                   // for fu
 ET*,                                                 // for u
 ET**,                                                // for D
 ET*>                                                 // for c
Program;

//typedef CGAL::Quadratic_program<real_t> Program;
typedef CGAL::Quadratic_program_solution<ET> Solution;

#define QPROG_EPS real_t(0.01)


Laplacian_DUAL2::Laplacian_DUAL2() : Laplacian_DUAL() {
}

Laplacian_DUAL2::~Laplacian_DUAL2() {
}

void Laplacian_DUAL2::clear(bool _c) {
    if (init) {
        cholmod_free_dense(&F_in, c);
        cholmod_free_dense(&F_ext, c);
        delete [] N0;

        // Call parent deletion
        Laplacian_DUAL::clear(_c);
    }

}

cholmod_sparse *Laplacian_DUAL2::constructLHS(Mesh *m) {
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

    // Calculate our differential coordinate matrices L
    Laplacian_DUAL::diffMatrix(m, L_trip);
    cholmod_sparse *L = cholmod_triplet_to_sparse(L_trip,L_trip->nnz,c);
    cholmod_free_triplet(&L_trip, c);       

    // Now compute the second order differential matrix. First we must compute a temporary first order
    // laplacian matrix Delta
    cholmod_dense *Delta = cholmod_allocate_dense(L->nrow, X->ncol, L->nrow, CHOLMOD_REAL, c);
    double alpha[2] = {1.0, 0.0};
    double beta[2] = {0.0, 0.0};
    cholmod_sdmult(L, 0, alpha, beta, X, Delta, c);

    // Now fill the triplet matrix with data
    cholmod_triplet *K_trip = cholmod_allocate_triplet(L->nrow, L->nrow, L->nrow*4, 0, CHOLMOD_REAL, c);
    secondDiffMatrix(K_trip, Delta);
    cholmod_sparse *K = cholmod_triplet_to_sparse(K_trip,K_trip->nnz,c);
    cholmod_free_triplet(&K_trip, c);
    cholmod_free_dense(&Delta, c);

    // Compute the product matrix KL
    cholmod_sparse *KL = cholmod_ssmult(K, L, 0, 1, 1, c);

    cholmod_triplet *KL_trip = cholmod_sparse_to_triplet(KL, c);    
    // Make sure we have enough space for the data in the triplet structure
    cholmod_reallocate_triplet(KL_trip->nnz+anchorList.size()+handleList.size(), KL_trip, c);
    //cholmod_print_triplet(KL_trip, "KL_trip", c);
    cholmod_free_sparse(&K, c);
    cholmod_free_sparse(&KL, c);
    cholmod_free_sparse(&L,c);

    int* KL_tripi = (int*) KL_trip->i;
    int* KL_tripj = (int*) KL_trip->j;
    double* KL_tripx = (double*) KL_trip->x;

    // Adding our constraints to the end of our matrix L
    list<int>::iterator it;
    int cnt;

    // Add our anchors first
    for (it = anchorList.begin(), cnt = dm.n_vertices(); it != anchorList.end(); ++it, ++cnt) {
        // Do a little dynamic resizing if necessary
        if (KL_trip->nnz+1 >= KL_trip->nzmax) {
            cholmod_reallocate_triplet(KL_trip->nnz + 64, KL_trip, c);
        }
        // Append our anchor
        KL_tripi[KL_trip->nnz] = cnt;
        KL_tripj[KL_trip->nnz] = *it;
        KL_tripx[KL_trip->nnz] = 1.0;
        KL_trip->nnz++;
    }

    // Now add our handles
    for (it = handleList.begin(); it != handleList.end(); ++it, ++cnt) {
        // Do a little dynamic resizing if necessary
        if (KL_trip->nnz+1 >= KL_trip->nzmax) {
            cholmod_reallocate_triplet(KL_trip->nnz + 64, KL_trip, c);
        }
        // Append our anchor
        KL_tripi[KL_trip->nnz] = cnt;
        KL_tripj[KL_trip->nnz] = *it;
        KL_tripx[KL_trip->nnz] = 1.0;
        KL_trip->nnz++;

        // Set our map value
        handleMap[*it] = cnt;
    }

    // Make L_trip into a sparse
    cholmod_sparse *ret_val = cholmod_triplet_to_sparse(
                    KL_trip, /* matrix to copy */
                    KL_trip->nnz,          /* allocate at least this much space in output matrix */
                    c);

    cholmod_free_triplet(&KL_trip, c);


    // Calculate the old normals of the dual faces (there is probably a better place to put this!)
    vector<double> N(3,0);
    double a;
    int i,j;
    N0 = new double[dm.n_faces()*3];
    DualMesh::FaceIter f_it;
    for (f_it = dm.faces_begin(), i=0; f_it != dm.faces_end(); ++f_it, ++i) {
        a = findNormal(f_it.handle(), &N);
        for (j=0; j<3; ++j) {
            N0[i+j*dm.n_faces()] = N[j];
        }
    }

    return ret_val;
}

/**
  * Construct the RHS of the problem. Note that this is different to the standard approach
  * as we have two force vectors Fin and Fext representing initial and external forces
  * respectively.
  */
void Laplacian_DUAL2::constructRHS(cholmod_sparse *L) {
    // Just allocate a big chunk for our dense matrix B, F_in and F_ext
    B = cholmod_allocate_dense(
            dm.n_vertices()+anchorList.size()+handleList.size(), // size_t nrow, /* # of rows of matrix */
            dim,    // size_t ncol, /* # of columns of matrix */
            dm.n_vertices()+anchorList.size()+handleList.size(), // size_t d, /* leading dimension */
            CHOLMOD_REAL,  // int xtype, /* CHOLMOD_REAL, _COMPLEX, or _ZOMPLEX */
            c);
    F_in = cholmod_allocate_dense(
            dm.n_vertices(), // size_t nrow, /* # of rows of matrix */
            dim,    // size_t ncol, /* # of columns of matrix */
            dm.n_vertices(), // size_t d, /* leading dimension */
            CHOLMOD_REAL,  // int xtype, /* CHOLMOD_REAL, _COMPLEX, or _ZOMPLEX */
            c);
    F_ext = cholmod_allocate_dense(
            dm.n_vertices(), // size_t nrow, /* # of rows of matrix */
            dim,    // size_t ncol, /* # of columns of matrix */
            dm.n_vertices(), // size_t d, /* leading dimension */
            CHOLMOD_REAL,  // int xtype, /* CHOLMOD_REAL, _COMPLEX, or _ZOMPLEX */
            c);

    double alpha[2] = {1.0, 0.0};
    double beta[2] = {0.0, 0.0};

    // Fill B with some initial values and constraints.
    cholmod_sdmult(L, 0, alpha, beta, X, B, c);

    // Copy some data from B into F_in for the initial forces of the model
    uint i;
    double *F_inx = (double*) F_in->x;
    double *Bx = (double*) B->x;
    for (i=0; i<F_in->ncol; ++i) {
        memcpy(&(F_inx[i*F_in->nrow]),&(Bx[i*B->nrow]), sizeof(double)*F_in->nrow);
    }

    // Set the data in F_ext to zero for now
    memset(F_ext->x, 0, sizeof(double)*F_ext->ncol*F_ext->nrow);

    // Set something in B
    updateB();
}

/**
  * Iteratively solve the second order dual Laplacian using
  */
void Laplacian_DUAL2::update() {
    // Sanity check
    if (!init) return;


    // Perform one loop, storing the current vertex positions
    cholmod_dense *_X = cholmod_copy_dense(X, c);

    updateDX();             // Update the dual vertices
    correct();              // Update the detail vectors according to local rotations
    updateB();             // Assemble the answer vector B
    Laplacian::update();    // Call the parent update function to solve for X

    double eps = epsilon(_X, X);
    int it = 0;
    int max_it = 100;
    double min_eps = 0.0001;
    fprintf(stderr,"\n - %d: epsilon(_X,X) = %f",it, eps);
    //char fname[256];
    double old_eps = eps;
    while ((eps > min_eps) && (it < max_it)  && (old_eps >= eps)) {
        old_eps = eps;
        cholmod_copy_dense2(X, _X, c);
        updateDX();         // Update the dual vertices
        correct();          // update the detail vectors according to local rotations
        updateB();             // Assemble the answer vector B
        Laplacian::update();    // Call the parent update function to solve for X
        eps = epsilon(_X, X);   // Compute the error
        ++it;
        fprintf(stderr,"\n - %d: epsilon(_X,X) = %f",it, eps);

        // For debuggering
        //sprintf(fname, "../data/debug_%d.obj", it);
        //dm.setX(DX);
        //dm.save(fname);
    }
}

void Laplacian_DUAL2::update_cuda() {
    // Sanity check
    if (!init) return;

    double p = 1.0;

    // Initialise cuda
    updateDX();
    CudaDual2::init_dual2_cu(DX->nrow, dm.n_faces(), B->nrow, dm.getTri(), (double*)DX->x, (double*)F_in->x, (double*)F_ext->x, N0, &(faceArea[0]));

    // Perform one loop, storing the current vertex positions
    cholmod_dense *_X = cholmod_copy_dense(X, c);

    updateDX();             // Update the dual vertices
    CudaDual2::correct_dual2_cu((double*)DX->x, (double*)B->x, p); // Update the detail vectors according to local rotations
    Laplacian::update();    // Call the parent update function to solve for X

    double eps = epsilon(_X, X);
    int it = 0;
    int max_it = 100;
    double min_eps = 0.0001;
    fprintf(stderr,"\n - %d: epsilon(_X,X) = %f",it, eps);
    //char fname[256];
    double old_eps = eps;

    while ((eps > min_eps) && (it < max_it)) { // && (old_eps >= eps)) {
        old_eps = eps;
        cholmod_copy_dense2(X, _X, c);
        updateDX();         // Update the dual vertices
        CudaDual2::correct_dual2_cu((double*)DX->x, (double*)B->x, p); // Update the detail vectors according to local rotations
        Laplacian::update();    // Call the parent update function to solve for X
        eps = epsilon(_X, X);   // Compute the error
        ++it;
        fprintf(stderr,"\n - %d: epsilon(_X,X) = %f",it, eps);

        // For debuggering
        //sprintf(fname, "../data/debug_%d.obj", it);
        //dm.setX(DX);
        //dm.save(fname);
    }
    // Delete GPU allocated memory
    CudaDual2::shutdown_dual2_cu();
}

/**
  * Compute the value of B by summing the forces.
  */
void Laplacian_DUAL2::updateB() {
    double *Bx = (double*)B->x;
    double *F_inx = (double*)F_in->x;
    double *F_extx = (double*)F_ext->x;
    uint i,j;
    for (i=0; i<F_in->nrow; ++i) {
        for (j=0; j<F_in->ncol; ++j) {
            Bx[i+j*B->nrow] = F_inx[i+j*F_in->nrow]+F_extx[i+j*F_in->nrow];
        }
    }
}

/**
  *
  */
void Laplacian_DUAL2::secondDiffMatrix(cholmod_triplet *K_trip, cholmod_dense *Delta) {
    // Our random access into the dual matrix
    int* Ki = (int*) K_trip->i;
    int* Kj = (int*) K_trip->j;
    real_t* Kx = (real_t*) K_trip->x;

    // Loop over all our dual facets
    Mesh::FaceIter f_it;
    Mesh::VertexIter v_it;
    Mesh::FaceVertexIter fv_it;
    int idx,i;
    vector<real_t> w_ij(3,real_t(0));

    for (f_it = dm.faces_begin(), v_it=dm.vertices_begin(), idx=0; f_it != dm.faces_end(); ++f_it, ++idx, ++v_it) {
        // Find the weights of the given simplex and its parent vertex
        findSecondWeights(Delta, v_it.handle(), f_it.handle(), &w_ij);

        // Assign the weights to our dual laplacian
        for (fv_it = dm.fv_iter(f_it.handle()),i=0;fv_it;++fv_it,++i) {
            Ki[K_trip->nnz] = idx;
            Kj[K_trip->nnz] = fv_it.handle().idx();
            Kx[K_trip->nnz] = w_ij[i];
            K_trip->nnz++;
        }

        // Let us hope that our barycentric coordinates sum to 1!
        Ki[K_trip->nnz] = idx;
        Kj[K_trip->nnz] = idx;
        Kx[K_trip->nnz] = -1.0;
        K_trip->nnz++;
    }
}

/**
 * Correct the matrix B by assigning to it the value -h_i * n_i
 */
void Laplacian_DUAL2::correct() {
    // Now we have the new vertices, we compute the normal for each face of
    // the dual complex, and update the matrix B with the result
    real_t *F_inx = (real_t*) F_in->x;
    real_t *F_extx = (real_t*) F_ext->x;
    int i,j;
    Mesh::FaceIter f_it;        // Face iterator
    vector<real_t> N(3,0);      // A normal compute from the current face
    vector<real_t> ON(3,0);     // Something in which to store the old normal
    vector<real_t> CP(3,0);     // Somewhere to store the cross product
    real_t d;                   // Somewhere to store the dot product
    real_t a;                   // The area of the current face
    real_t scale;               // The scaling factor based on face area
    vector<real_t> F_0(3,0);    // A vector storing F_in before rotation
    vector<real_t> F(3,0);      // A vector storing F after rotation

    cholmod_dense *R = cholmod_allocate_dense(3,3,3,CHOLMOD_REAL,c);
    real_t *Rx = (real_t*) R->x;
    real_t phi;

    for (f_it = dm.faces_begin(), i=0; f_it != dm.faces_end(); ++f_it, ++i) {
        // Compute the rotation matrix from the input normal
        a = findNormal(f_it.handle(), &N);
        for (j=0; j<3; ++j) {
            F_0[j] = F_inx[i+j*F_in->nrow];
            ON[j] = N0[i+j*dm.n_faces()];
        }
        DOT3(d,N,ON); CROSS3(CP,N,ON); NORMALIZE3(CP);
        //fprintf(stderr, "\nN=[%f,%f,%f], ON=[%f,%f,%f], CP=[%f,%f,%f]",
        //        N[0],N[1],N[2],ON[0],ON[1],ON[2],CP[0],CP[1],CP[2]);

        phi = acos(d);
        if (phi > 0.00001) {
            //fprintf(stderr, "\nphi=%f, N X N0 = [%f,%f,%f]",phi,CP[0],CP[1],CP[2]);
            rotMat(R, &CP, phi);
            //cholmod_print_dense(R, "Rotation", c);
            // Apply the rotation to the F_in
            MATVECMULT3(F, Rx, F_0);
        } else {
            for (j=0; j<3; ++j) F[j] = F_0[j];
        }

        scale = sqrt(a/faceArea[i]);

        // Copy the result back into the F_in and F_ext
        for (j = 0; j < 3; ++j) {
            F_inx[i + j*F_in->nrow ] = scale * F[j];
            F_extx[i + j*F_ext->nrow ] = scale*F_extx[i+j*F_ext->nrow];

            // Update the normal for this face
            N0[i+j*dm.n_faces()] = N[j];
        }
        // Update the face area
        faceArea[i] = a;
    }
    // Delete allocated memory
    cholmod_free_dense(&R, c);
}

/**
  * Do something clever to the anchors here.
  */
void Laplacian_DUAL2::transformHandles(cholmod_dense *M) {
    Laplacian::transformHandles(M);
    return;
}


/**
  * Specify an external force on the surface by applying a radial decay function
  * and computing the forces to all the points on the mesh. Uses a gaussian function
  * to compute the value to scale the force at each vertex.
  */
void Laplacian_DUAL2::applyRadialForce(const double *pos,
                                        const double *force,
                                        const double radius) {


    //fprintf(stderr, "\nLaplacian_DUAL2::applyRadialForce([%f,%f,%f], [%f,%f,%f], %f)",
    //        pos[0],pos[1],pos[2],force[0],force[1],force[2],radius);

    double inv_rad = 1.0 / (radius*radius);
    double *F_x = (double*)F_ext->x;
    uint i,j;
    double dist, wt;

    Mesh::VertexIter vit;
    for (vit = dm.vertices_begin(), i=0; vit != dm.vertices_end(); ++vit, ++i) {
        // Compute the sqrd distance from pos to the current point
        dist = 0.0;
        for (j=0; j<3; ++j)
            dist += (dm.point(vit.handle())[j] - pos[j])*(dm.point(vit.handle())[j] - pos[j]);


        // Update the external force
        wt = exp(-dist*inv_rad);

        //fprintf(stderr, "\n - wt(%d) = %f, pt=[%f,%f,%f], dist=%f, inv_rad=%f", i, wt,
        //        dm.point(vit.handle())[0],dm.point(vit.handle())[1],dm.point(vit.handle())[2],
        //        dist, inv_rad);

        for (j=0; j<3; ++j)
            F_x[i+j*F_ext->nrow] += wt * force[j]; // / faceArea[i];
    }
}

/**
  * Using CGAL for quadratic optimization. This should be as follows:
  * CGAL solves problem:
  * min x'Dx + c'x + c_0
  * s.t. Ax >=< b, l<=x<=u
  *
  * P = [y_{j,1}, y_{j,2}, y_{j,3}]
  * k = [|y_i - y_{j,1}|; |y_i - y_{j,2}|; |y_i - y_{j,3}|]
  * a is some weight (something like 0.01)
  * D = P'P + ak'k
  * c = (-2y_i'A)'
  * c_0 = y_i'y_i
  * Also we add an equality constraint to ensure that the weights sum to one (A_1=[1,1,1], b=[1])
  */
void Laplacian_DUAL2::findSecondWeights(cholmod_dense *B, Mesh::VertexHandle vh, Mesh::FaceHandle fh, vector<real_t> *wt) {    
    // The indices of the face handle must be retrieved
    int i,j;
    int f_i = vh.idx();
    vector<int> f_ij(3,0);
    Mesh::FaceVertexIter fvit;
    for (fvit = dm.fv_iter(fh),i=0; fvit; ++fvit,++i) {
        f_ij[i] = fvit.handle().idx();
    }

    // First build up matrix P and vector k, which we'll use for D and c
    double P[9];
    double D[9];
    ET A_col[] = {ET(1)};
    ET* A[] = {A_col, A_col, A_col};
    double k[] = {double(0),double(0),double(0)};
    ET q[] = {ET(0),ET(0),ET(0)};

    bool fl[] = {true,true,true};
    bool fu[] = {true,true,true};
    Relation r(CGAL::EQUAL);
    ET l[] = {QPROG_EPS,QPROG_EPS,QPROG_EPS};
    ET u[] = {ET(1),ET(1),ET(1)};
    ET b[] = {ET(1)};

    double *Bx = (real_t*) B->x;
    for (i=0; i<3; ++i) {        
        for (j=0; j<3; ++j) {
            // P has the points as columns
            P[j+i*3] = Bx[f_ij[i] + j*B->nrow];
            k[i] += (Bx[f_i + j*B->nrow] - P[j+i*3]) * (Bx[f_i + j*B->nrow] - P[j+i*3]);
            q[i] += ET(Bx[f_i + j*B->nrow]*P[j+i*3]);
        }
        k[i] = sqrt(k[i]);
        q[i] *= -ET(2);
    }

    // Solve for kkt manually!
    D[0] = k[0]*k[0];
    D[4] = k[1]*k[1];
    D[8] = k[2]*k[2];
    D[1] = D[3] = k[0]*k[1];
    D[2] = D[6] = k[0]*k[2];
    D[5] = D[7] = k[1]*k[2];

    // Test the values of D are not all close to zero, which can arise if the thing is flat
    bool is_feasible = false;
    for (i=0; i<9; ++i) {
        is_feasible = is_feasible || (D[i] > 1.0e-10);
    }
    if (!is_feasible) {
        fprintf(stderr, "\nLaplacian_DUAL2::findSecondWeights() - D = zeros(3) - region flat!");
        for (i=0; i<3; ++i) wt->at(i) = 1.0 / 3.0;
        return;
    }

    //fprintf(stderr, "\nD=[%f,%f,%f;%f,%f,%f;%f,%f,%f]",D[0],D[1],D[2],D[3],D[4],D[5],D[6],D[7],D[8]);

    // Solve for D = PtP + akk'
    int M = 3; // Rows of A
    int N = 3; // Cols of B
    int K = 3; // Cols of A and Rows of B
    int LDA = M; // The leading dimension of matrix A
    int LDB = K; // Leading dimension of matrix B
    int LDC = M; // Leading dimension of matrix C
    double alpha = 1.0; // A is scaled by this much
    double beta = QPROG_EPS; // C is scaled by this much (change this parameter to change the weight of edges!)
    char TRANSA = 'T'; // Transpose A
    char TRANSB = 'N'; // Do not transpose B
    dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha, P, &LDA, P, &LDB, &beta, D, &LDC);

    // Copy the result into _D

    ET _D_row1[] = {ET(D[0])};
    ET _D_row2[] = {ET(D[1]), ET(D[4])};
    ET _D_row3[] = {ET(D[2]), ET(D[5]), ET(D[8])};
    ET* _D[] = {_D_row1, _D_row2, _D_row3};

    Program qp(3, 1, A, b, r, fl, l, fu, u, _D, q);

    // now solve the problem
    //dumpQP(qp);
    Solution s = CGAL::solve_quadratic_program(qp, ET());
    if (!s.is_void() && s.is_optimal()) {
        Solution::Variable_value_iterator vit = s.variable_values_begin();
        for (i=0;vit != s.variable_values_end(); ++vit,++i) {
            wt->at(i) = to_double(*vit);
        }
        //cerr << "\nwt=["<<wt->at(0)<<","<<wt->at(1)<<","<<wt->at(2)<<"]";
    } else {
        //cgal_failed!        
        cerr << "\nCGAL failed! - ";
        if (s.is_infeasible()) cerr << "INFEASIBLE! ";
        if (s.is_unbounded()) cerr << "UNBOUNDED!";        
        for (i=0; i<3; ++i) wt->at(i) = 1.0 / 3.0;
    }
}

/** ROTMAT Construct the rotation matrix out of an up vector and an angle of
  * rotation using eulers method
  * http://mathworld.wolfram.com/EulerParameters.html.
  */
void Laplacian_DUAL2::rotMat(cholmod_dense *R, vector<real_t> *N, real_t phi) {
    real_t sin_phi = sin(0.5*phi);
    real_t e[] = {cos(0.5*phi), N->at(0)*sin_phi, N->at(1)*sin_phi, N->at(2)*sin_phi};
    real_t ee[] = {e[0]*e[0], e[1]*e[1], e[2]*e[2], e[3]*e[3]};

    // Create our return type
    real_t *Rx = (real_t*) R->x;

    Rx[0] = ee[0] + ee[1] - ee[2] - ee[3];
    Rx[1] = real_t(2)*(e[1]*e[2] - e[0]*e[3]);
    Rx[2] = real_t(2)*(e[1]*e[3] + e[0]*e[2]);
    Rx[3] = real_t(2)*(e[1]*e[2] + e[0]*e[3]);
    Rx[4] = ee[0] - ee[1] + ee[2] - ee[3];
    Rx[5] = real_t(2)*(e[2]*e[3] - e[0]*e[1]);
    Rx[6] = real_t(2)*(e[1]*e[3] - e[0]*e[2]);
    Rx[7] = real_t(2)*(e[2]*e[3] + e[0]*e[1]);
    Rx[8] = ee[0] - ee[1] - ee[2] + ee[3];
}

