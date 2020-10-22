#ifndef LAPLACIAN_DUAL2_H
#define LAPLACIAN_DUAL2_H

#include "laplacian_DUAL.h"
#include "laplacian_DUAL2_cuda.cuh"

/**
  * An implementation of the second order dual laplacian operator.
  */
class Laplacian_DUAL2 : public Laplacian_DUAL {
public:
    Laplacian_DUAL2();
    ~Laplacian_DUAL2();

    /// Compute the diffMatrix
    void diffMatrix(cholmod_triplet */*A*/);

    /// Solve the system using a special approach
    virtual void update();
    virtual void update_cuda();

    /// Special handling of anchor transforms
    virtual void transformHandles(cholmod_dense */*M*/);

    virtual void applyRadialForce(const double *pos,
                                  const double *force,
                                  const double radius);

protected:
    virtual cholmod_sparse *constructLHS(Mesh *);       //< Construct the laplacian matrix L
    virtual void constructRHS(cholmod_sparse *);        //< Construct B

    /// Compute the second order laplacian weights for the matrix K
    void secondDiffMatrix(cholmod_triplet */*K_trip*/, cholmod_dense */*Delta*/);
    void findSecondWeights(cholmod_dense */*B*/, Mesh::VertexHandle /*vh*/, Mesh::FaceHandle /*fh*/, vector<double> */*wt*/);

    /// Compute the answer vector B by summing forces
    void updateB();

    /// Correct the force vectors according to the dual mesh deformation
    virtual void correct();

    /// Delete allocated memory
    virtual void clear(bool);
private:
    /// Initial and external force vectors
    cholmod_dense *F_in, *F_ext;

    /// Store the initial normals
    double *N0;

    /// Rotation matrix computation, stored in R
    void rotMat(cholmod_dense */*R*/, vector<double> */*N*/, double /*phi*/);
};

#endif //LAPLACIAN_DUAL2_H
