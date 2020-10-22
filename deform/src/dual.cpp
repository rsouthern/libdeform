#include "dual.h"
#include <vector>
#include <map>

using namespace std;

/**
  * Construct a dual mesh from the input mesh class.
  */
DualMesh::DualMesh(Mesh *m) : Mesh() {
    // Initialize cholmod
    c = CholmodSharedSingleton::Instance()->getCholmodCommon();
    if (m != NULL) construct(m);
}

/**
  * Destroy this object
  */
DualMesh::~DualMesh() {
    cholmod_free_sparse(&D, c);
    CholmodSharedSingleton::stop();
}

/**
  * Call this to construct the mesh (need to work out how to clean this!)
  */
void DualMesh::construct(Mesh *m) {
    // Create a vertex handle for each of the vertices so we can link them together
    vector<VertexHandle> vhandles(m->n_faces());

    ConstFaceIter f_it;

    // Create our cholmod triplet to prepare our sparse data structure
    cholmod_triplet *Dt = cholmod_allocate_triplet(m->n_faces(), m->n_vertices(), m->n_faces()*3, 0, CHOLMOD_REAL, c);
    int *Di = (int*) Dt->i;
    int *Dj = (int*) Dt->j;
    real_t *Dx = (real_t*) Dt->x;
    real_t dval = 1.0 / 3.0;

    // First create all our dual vertices
    for (f_it=m->faces_begin(); f_it!=m->faces_end(); ++f_it) {
        // First create the entry in the dual matrix
        Point pt(0.0,0.0,0.0);
        for (ConstFaceVertexIter cfv_it = m->cfv_iter(f_it.handle()); cfv_it; ++cfv_it) {
            Di[Dt->nnz] = f_it.handle().idx();
            Dj[Dt->nnz] = cfv_it.handle().idx();
            Dx[Dt->nnz] = dval;
            Dt->nnz++;
            pt += m->point(cfv_it);
        }
        // Now create the dual vertex in our new mesh
        vhandles[f_it.handle().idx()] = add_vertex(Point(dval*pt[0], dval*pt[1], dval*pt[2]));
    }

    // A structure to store the vertex handles of our dual face
    vector<Mesh::VertexHandle> face_vhandles;

    // Now create all our dual vertices using face iterators
    for (Mesh::ConstFaceIter f_it=m->faces_begin(); f_it!=m->faces_end(); ++f_it) {
        face_vhandles.clear();
        // Iterate over the half edges of the face
        for (Mesh::ConstFaceFaceIter ff_iter = m->cff_iter(f_it); ff_iter; ++ff_iter) {
            // Determine if there is a face on the other side
            face_vhandles.push_back(vhandles[ff_iter.handle().idx()]);

        }
        add_face(face_vhandles);
    }
    D = cholmod_triplet_to_sparse(Dt, Dt->nnz, c);
    cholmod_free_triplet(&Dt, c);

    // Update the normals and the storage vector X
    update_normals();
    updateX();
    updateTri();
}


/**
  * Return a pointer to our Dual Operator matrix.
  */
cholmod_sparse *DualMesh::getD() {
    return D;
}

