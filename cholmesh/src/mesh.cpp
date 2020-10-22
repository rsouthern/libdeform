#include "mesh.h"
#include <cstdlib>

/**
  * Initialise the mesh class. Tells cholmod we're using it.
  */
Mesh::Mesh() : MeshBase() {
    // Add vertex normals as default property (ref. previous tutorial)
    request_vertex_normals();

    // Add face normals as default property
    request_face_normals();

    // Get a cholmod workspace
    c = CholmodSharedSingleton::Instance()->getCholmodCommon();

    // Initialise the data to NULL
    X = NULL;
    Tri = NULL;
}

/**
  * Destructor. Tells cholmod we're no longer using it.
  */
Mesh::~Mesh() {
    if (X != NULL) cholmod_free_dense(&X, c);
    if (Tri != NULL) delete [] Tri;
    CholmodSharedSingleton::stop();
}

/**
  * A wrapper for the OpenMesh IO function. Also initialises the normals.
  */
int Mesh::load(const char* fname) {
    int res = OpenMesh::IO::read_mesh(*this, fname);
    if (res) update_normals();
    updateX();
    updateTri();
    // A little test to see that the normals were working - they are! ;)
    /*
    // move all vertices one unit length along it's normal direction
    for (Mesh::VertexIter v_it = vertices_begin(); v_it != vertices_end(); ++v_it) {
      set_point( v_it, point(v_it)+ ((float) rand())/((float) RAND_MAX) * normal(v_it) );
    }*/
    return res;
}

/**
  * Save this mesh to an output file using OpenMesh IO
  */
int Mesh::save(const char* fname) {
    return OpenMesh::IO::write_mesh(*this, fname);
}

/**
  * Return the local copy of the data X. This class owns the member X.
  */
cholmod_dense *Mesh::getX() {
    return X;
}

/**
  * Return the triangulation data structure
  */
uint *Mesh::getTri() {
    return Tri;
}

/**
  * Convert the point data to a dense matrix for CHOLMOD (data is not managed by this class)
  */
void Mesh::updateX() {
    if (X != NULL) cholmod_free_dense(&X, c);

    // Allocate our memory
    X = cholmod_allocate_dense(n_vertices(), 3, n_vertices(), CHOLMOD_REAL, c);
    if (c->status != 0) {
        fprintf(stderr, "\nMesh::toDense() - Could not allocate X[%d,%d], status=%d!", n_vertices(), 3, c->status);
        return;
    }

    // Iterate through all vertices and retrieve the point data
    int i,j;
    Mesh::VertexIter v_it;
    for (v_it = vertices_begin(), i=0; v_it != vertices_end(); ++v_it, ++i) {
        for (j=0;j<3;++j) {
            ((real_t*) X->x)[i + j*n_vertices()] = point(v_it)[j];
        }
    }
}

/**
  * The Tri vector contains the vertex data for the triangle.
  */
void Mesh::updateTri() {
    if (Tri != NULL) delete [] Tri;

    // Allocate memory
    Tri = new uint [n_faces()*3];

    // Iterate over the triangles and assign the indices to the data structure
    Mesh::FaceIter f_it;
    Mesh::FaceVertexIter fv_it;
    int i,j;
    int N = n_faces();
    for (f_it = faces_begin(),i=0; f_it != faces_end(); ++f_it,++i) {
        for (fv_it = fv_iter(f_it.handle()),j=0; fv_it; ++fv_it,++j) {
            Tri[i + j*N] = fv_it.handle().idx();
        }
    }
}

/**
  * Set our point data from a dense matrix structure.
  */
void Mesh::setX(cholmod_dense *_X) {
    // Sanity check
    if ((_X->nrow != n_vertices()) || (_X->ncol != 3)) {
        fprintf(stderr, "\nMesh::setX(X) - Data structure is incorrect size! Expected [%d,%d], got [%d,%d]",
                (int) n_vertices(), 3, (int) _X->nrow, (int) _X->ncol);
        return;
    }

    // Set our local data structure
    if (X != NULL) cholmod_free_dense(&X, c);
    X = cholmod_copy_dense(_X, c);
    real_t *Xx = (real_t*) X->x;

    // Now update the mesh structure positions
    Mesh::VertexIter v_it;
    int j;
    for (v_it = vertices_begin(); v_it != vertices_end(); ++v_it) {
        for (j=0; j<3; ++j)
            point(v_it)[j] = Xx[v_it.handle().idx() + j*n_vertices()];
    }
}

