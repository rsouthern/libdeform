#ifndef _GS_VECMATH_H
#define _GS_VECMATH_H

#include "gs_Defines.h"
#include <math.h>

/** Define my own simple math routines. Assuming all 
  * vectors are of type real[3].
  */
namespace VecMath {
    /// Determine the squared distance between two 
    real sqrDist(real* p1, real* p2);
    
    /// Determine the actual distance (sqrt'ed)
    real dist(real* p1, real* p2);
    
    /// Normalise a vector
    void normalise(real* v);
    
    /// Add two 3-vectors
    void add(real* ans, real* v1, real* v2);
    
    /// Subtract two 3-vectors
    void sub(real* ans, real* v1, real* v2);
    
    /// Determine the cross product of 2 3-vectors
    void cross(real* ans, real* v1, real* v2);
    
    /// Return the dot product of 2 3-vectors
    real dot(real* v1, real* v2);
    
    /// Determine the l2 norm of a 3 vector
    real l2(real*);
    
    /// Computes the area of given triangle given 3 pts
    real triArea(real* p1, real* p2, real* p3);
    
};

#endif //_GS_VECMATH_H

