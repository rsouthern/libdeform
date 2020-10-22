#ifndef _LAPLACIAN_DUAL_H
#define _LAPLACIAN_DUAL_H

#include "laplacian.h"
#include "dual.h"
#include <vector>
#include "laplacian_DUAL_cuda.cuh"

using namespace std;

typedef Mesh::VertexHandle VHandle;
typedef Mesh::FaceHandle FHandle;

// Cross product
#define CROSS3(dest,v1,v2) \
          dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
          dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
          dest[2]=v1[0]*v2[1]-v1[1]*v2[0];
// Dot product
#define DOT3(dest,v1,v2) \
          dest=v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
// Assume column major ordering!
#define MATVECMULT3(dest,M,v) \
          dest[0]=M[0]*v[0]+M[3]*v[1]+M[6]*v[2];\
          dest[1]=M[1]*v[0]+M[4]*v[1]+M[7]*v[2];\
          dest[2]=M[2]*v[0]+M[5]*v[1]+M[8]*v[2];
// Subtract two vectors
#define SUB3(dest, v1, v2) \
          dest[0]=v1[0]-v2[0];\
          dest[1]=v1[1]-v2[1];\
          dest[2]=v1[2]-v2[2];
#define NORMALIZE3(v) \
          double inv_norm=1.0 / sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); \
          v[0] *= inv_norm; \
          v[1] *= inv_norm; \
          v[2] *= inv_norm;

class Laplacian_DUAL : public Laplacian {
public:
    /// Constructor
    Laplacian_DUAL();

    /// Destructor
    virtual ~Laplacian_DUAL();

    /// Build the laplacian matrix and solve for the inverse
    virtual void initialise(Mesh */*m*/);   

    /// Restore the state of the geometry iteratively
    virtual void update();

    /// Throw caution to the wind and update the detail vectors on the GPU!
    virtual void update_cuda();

    /// Overload this function so it will get called
    virtual bool needsRestore() {return true;}

protected:
    /// vector used to store magnitude of local differential feature
    vector<double> correction;

    /// vector used to store face area of individual feature
    vector<double> faceArea;

    /// The dual mesh which we'll keep a copy of to update surface normals
    DualMesh dm;

    /// Store our dual vertices in DX - computed each time the mesh is modified
    cholmod_dense *DX;        

    /// Inherited from parent
    void clear(bool /*clear constraints*/ = false);

    /// Find the weights of the dual graph laplacian
    void findWeights(VHandle, FHandle, vector<double> */*w_ij*/, double */*h_i*/, double */*area_i*/);

    /// Finds the maximum distance between two point sets
    double epsilon(cholmod_dense */*first*/, cholmod_dense */*second*/);

    /// Fix B with the updated normals and the correction value
    virtual void correct();

    /// Find the normal to the face and set the parameter. Also return face area.
    double findNormal(Mesh::FaceHandle /*fh*/, vector<double> */*normal*/);

    /// Update the dual vertices
    void updateDX();

    // The following functions are called during initialisation
    virtual cholmod_sparse *constructLHS(Mesh *);       //< Construct the laplacian matrix L
    virtual void constructRHS(cholmod_sparse *);        //< Construct B
    virtual void diffMatrix(Mesh */*m*/, cholmod_triplet */*A*/); //< Construct a matrix of differential coordinates
};

#endif //_LAPLACIAN_DUAL_H


