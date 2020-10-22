#include "laplacian_COT.h"
#include <vector>

using namespace std;

Laplacian_COT::Laplacian_COT() : Laplacian() {
}

Laplacian_COT::~Laplacian_COT() {
}

/**
 * Construct a matrix of differential coordinates for this complex. Note that
 * it doesn't have to be a manifold - it could be a full space complex.
 * The form of our system is Ax=B where A is a n*n matrix and B is a n*m matrix,
 * where n is the number of vertices and m is the dimension of the vertices.
 */
void Laplacian_COT::diffMatrix(Mesh *m, cholmod_triplet *A) {
    if (init) return;

    // Finding the differential coordinates for each vertex
    int idx;
    vector<double> w_ij;
    double sum_w_ij;

    int* Ai = (int*) A->i;
    int* Aj = (int*) A->j;
    double* Ax = (double*) A->x;

    uint oRsize;
    Mesh::VertexIter v_it;
    Mesh::VertexVertexIter vv_it;
    Mesh::VertexOHalfedgeIter voh_it;

    for (v_it = m->vertices_begin(); v_it != m->vertices_end(); ++v_it) {
        int v_idx = v_it.handle().idx();

        // Count the size of the one ring (is there a better way to do this?
        for (vv_it = m->vv_iter(v_it.handle()), oRsize=0; vv_it; ++vv_it, ++oRsize);

        // Allocate some space for our local weights (a ring size of zero is bad - skip it)
        if (oRsize > 0) {
            w_ij.resize(oRsize);
        } else {
            fprintf(stderr, "\nLaplacian_COT::diffMatrix() - Vertex %d has no neighbours!", v_idx);
            continue;
        }

        sum_w_ij = 0.0;

        // Iterate over all one ring half edges to compute the diff weights
        for (voh_it = m->voh_iter(v_it.handle()), idx = 0; voh_it; ++voh_it, ++idx) {
            // Compute the weight and add it to the total
            w_ij[idx] = diffWeight(m, voh_it.handle());
            sum_w_ij += w_ij[idx];
        }	

        double inv_w = 1.0 / sum_w_ij;
        //fprintf(stderr, "\nLaplacian_COT::diffMatrix() - inv_w = %f", inv_w);

        // For each vertex in the one ring
        for (vv_it = m->vv_iter(v_it.handle()), idx=0; vv_it; ++vv_it,++idx) {
            // Storing the fraction into the matrix A
            int vv_idx = vv_it.handle().idx();
            Ai[A->nnz] = v_idx;
            Aj[A->nnz] = vv_idx;
            Ax[A->nnz] = w_ij[idx]*inv_w;
            A->nnz++;
        }
        // Add our i entry in this row
        Ai[A->nnz] = v_idx;
        Aj[A->nnz] = v_idx;
        Ax[A->nnz] = -1.0;
        A->nnz++;
    }
}

/**
  * Returns the cotangent weights, that is
  * w_ij = 1/2 ( cot(a) + cot(b) )
  * where i and j are vertex indices on the edge eminating from vertex i and a and b
  * are the angles opposite to that edge.
  * We use the Half edge structure to navigate the mesh - which makes this computation trivially easy!
  */
double Laplacian_COT::diffWeight(Mesh *m, HEHandle he) {
    double a = m->calc_sector_angle(m->next_halfedge_handle(he));
    double b = m->calc_sector_angle(m->next_halfedge_handle(m->opposite_halfedge_handle(he)));
    double ret_val = 0.5 * (tan ( M_PI * 0.5 - a) + tan ( M_PI * 0.5 - b));
    //fprintf(stderr, "\nLaplacian_COT::diffWeight() - a=%f, b=%f, returning %f", a, b, ret_val);
    return ret_val;
}

