#ifndef MESH_H
#define MESH_H

#include <stdio.h>
#include <string.h>
#include <list>
#include <vector>
#include <algorithm>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <defines.h>
#include "cholmod_shared.h"

// The triangle mesh class we use in this project
typedef OpenMesh::TriMesh_ArrayKernelT<> MeshBase;

using namespace std;


/**
  * A mesh class containing face and vertex data
  */
class Mesh : public MeshBase {
public:
    Mesh();
    virtual ~Mesh();

    virtual int load(const char*);      //< Load the mesh from disk
    virtual int save(const char*);      //< Save the mesh to disk

    void setX(cholmod_dense */*X*/);    //< Set our point data from external data

    cholmod_dense *getX();              //< Retrieve the dense matrix X
    uint *getTri();                     //< Retrieve the triangulation data in the form of a matrix

    void updateTri();                   //< Update the triangle data structure Tri
    void updateX();                     //< Update the dense matrix X
protected:
    cholmod_common *c;                  //< Store the common workspace here
    cholmod_dense *X;                   //< Stored column ordering data for points
    uint *Tri;                           //< Stored triangle indices
};

#endif // MESH_H
