#include "cusp_tools.cuh"
/*
extern "C" bool cholmod_factor_to_cusp(cholmod_factor *fac, dMat *dM) {

    // We need access to the CHOLMOD workspace
    cholmod_common* c = CholmodSharedSingleton::Instance()->getCholmodCommon();

    // First convert the factor to a sparse matrix and transpose it
    cholmod_sparse A;

    // - Copy the factor into a manually created sparse structure. Note that this is a shallow copy,
    //   so no memory needs to be erased after this function.
    A.nrow = fac->n; A.ncol = fac->n; A.nzmax = fac->nzmax; A.p = fac->p; A.i = fac->i;
    A.nz = fac->nz; A.x = fac->x; A.z = fac->z; A.stype = 0; A.itype = fac->itype;
    A.xtype = fac->xtype; A.dtype = fac->dtype; A.sorted = TRUE; A.packed = false;

    cholmod_print_factor(fac, "FAC", c);
    cholmod_print_sparse(&A, "A from a factor", c);

    // - Make the transpose matrix packed! The CSR format needs to be packed.
    cholmod_sparse *A_trans = cholmod_allocate_sparse(fac->n, fac->n, fac->nzmax, 0, 1, 1, CHOLMOD_REAL, c);

    // - Perform the unsymmetric transpose. Hopefully it will be packed (because A_trans->packed == 1).
    cholmod_transpose_unsym(&A, 1, NULL, NULL, 0, A_trans, c);



    // Now we have transposed the factor, we can write it into the output structure
    hMat hM(A_trans->ncol ,A_trans->nrow ,A_trans->nzmax );

    // - Copy the column pointers into the row offsets
    memcpy(&(hM.row_offsets[0]), (int*) A_trans->p, (A_trans->nrow+1)*sizeof(int));

    // - Copy the row indices into the column indices
    memcpy(&(hM.column_indices[0]), (int*) A_trans->i, (A_trans->nzmax)*sizeof(int));

    // - Copy the values
    memcpy(&(hM.values[0]), (real_t*) A_trans->x, sizeof(real_t)*A_trans->nzmax);

    // Create a return value and return it (the CUSP equality operator knows to copy the data apparently)
    (*dM) = hM;

    // Delete allocated memory
    cholmod_free_sparse(&A_trans, c);


    // Return success!
    return true;
}*/
