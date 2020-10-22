#ifndef DUAL_H
#define DUAL_H

#include <list>
#include <vector>
#include <algorithm>
#include "defines.h"
#include "mesh.h"

using namespace std;

/**
  * A class containing the matrices necessary to maintain a dual mesh operator
  */
class DualMesh : public Mesh {
public:
    DualMesh(Mesh *m = NULL);   //< Only allowable constructor
    virtual ~DualMesh();        //< Destructor
    cholmod_sparse *getD();     //< Retrieve the dual operator matrix
    virtual void construct(Mesh *m); //< Construct this mesh from a primal mesh

private:
    cholmod_common *c;          //< Cholmod workspace
    cholmod_sparse *D;          //< The sparse dual operator matrix (size nF*nV, type float)
};

#endif // DUAL_H
