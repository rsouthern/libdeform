#ifndef _LAPLACIAN_COT_H
#define _LAPLACIAN_COT_H

#include "laplacian.h"

typedef Mesh::HalfedgeHandle HEHandle;

/**
 * Cotangent weights for Laplace Beltrami operator.
 */
class Laplacian_COT : public Laplacian {
public:
    /// Constructor
    Laplacian_COT();

    /// Destructor
    virtual ~Laplacian_COT();

protected:
    /// Construct a matrix of differential coordinates from the given complex
    virtual void diffMatrix(Mesh *m, cholmod_triplet */*A*/);

    /// Compute the local weight between two vertices by using the outgoing halfedge from v_i to v_j
    virtual double diffWeight(Mesh */*m*/, HEHandle /*he*/);
};

#endif //_LAPLACIAN_COT_H


