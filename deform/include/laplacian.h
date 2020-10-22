#ifndef _LAPLACIAN_H
#define _LAPLACIAN_H

#include <stdio.h>
#include <time.h>
#include <list>
#include <algorithm>
#include "cholmod_shared.h"
#include "mattools.h"
#include <mesh.h>
#include "defines.h"


using namespace std;

class Laplacian {
public:
    /// Default constructor and destructor
    Laplacian();
    virtual ~Laplacian();

    /// Clean up our laplacian to be used again
    virtual void clear(bool /*clear constraints*/ = false);

    /// Add/remove a constraint or anchor
    void addHandle(int i);
    void addAnchor(int i);

    /// Main interface - retrieve a pointer to deformed points
    cholmod_dense *getX();

    /// Build the laplacian matrix and solve for the inverse
    virtual void initialise(Mesh */*m*/);

    /// Deform the current laplacian based on the transformation of handles
    virtual void transformHandles(cholmod_dense */*transformation matrix*/ = NULL);
    
    /// Apply some post-process (optimisation) after the edit
    virtual void restore() {}

    /// Override this function if some postprocess is required
    virtual bool needsRestore() {return false;}

    /// Compute the answer vector for X
    virtual void update();

    /// Leave the option available for hard core programming!
    virtual void update_cuda() {
        update();
    }

    /// Dump the contents of this structure to screen
    virtual void print();

    /// Save the constraints to a file
    void saveConstraints(FILE */*fid*/);

    /// Load the constraints from a file
    void loadConstraints(FILE */*fid*/);

protected:
    /// The dimension of the vertices
    int dim;

    /// The number of vertices
    int vertNum;

    /// A variable storing whether this has been initialised
    bool init;

    /// A list of anchors and constraints
    list<int> anchorList, handleList;

    /// A map from handles to the row index L
    map<int,int> handleMap;

    /// Our answer vector which is used when we are deforming
    cholmod_dense *X;

    /// The matrix of constraint positions and differential coordinates
    cholmod_dense *B;

    /// The matrix which B must be multiplied by before solving for X
    cholmod_sparse *L_trans;

    /// A preallocated matrix for L'B -- used in deform()
    cholmod_dense *L_trans_B;

    /// The factorisation which only has to be done once the deformation begins
    cholmod_factor *fac;

    /// The workspace in which our cholmod objects live
    cholmod_common *c;

    // The following functions are called during initialisation
    virtual cholmod_sparse *constructLHS(Mesh *);       //< Construct the laplacian matrix L
    virtual void constructRHS(cholmod_sparse *);        //< Construct B
    virtual void prepareProblem(cholmod_sparse *);      //< Factorise the problem
    virtual void diffMatrix(Mesh */*m*/, cholmod_triplet */*A*/); //< Construct a matrix of differential coordinates

    /// For debugging
    void print_mem(const char* s) {printf("%s, %d / %d",s, (int) c->memory_inuse, (int) c->memory_usage);}

    /// Test the LHS matrix to evaluate the singular values
    virtual bool isMatValid(cholmod_sparse *);
};

#endif //_LAPLACIAN_H
