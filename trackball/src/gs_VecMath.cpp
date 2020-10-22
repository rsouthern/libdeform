#include "gs_VecMath.h"

/// Determine the squared distance between two 
real VecMath::sqrDist(real* p1, real* p2) {
    return  (p2[0]-p1[0])*(p2[0]-p1[0]) + 
	    (p2[1]-p1[1])*(p2[1]-p1[1]) + 
	    (p2[2]-p1[2])*(p2[2]-p1[2]);
}
    
/// Determine the actual distance (sqrt'ed)
real VecMath::dist(real* p1, real* p2) {
    return sqrt(sqrDist(p1, p2));
}    
    
/// Normalise a vector
void VecMath::normalise(real* v) {
    real length = l2(v);
    if (length != 0.0) {
	v[0]/=length; v[1]/=length; v[2]/=length;
    }	
}

    
/// Add two 3-vectors
void VecMath::add(real* ans, real* v1, real* v2) {
    ans[0] = v2[0] + v1[0];
    ans[1] = v2[1] + v1[1];
    ans[2] = v2[2] + v1[2];
}


/// Subtract two 3-vectors
void VecMath::sub(real* ans, real* v1, real* v2) {
    ans[0] = v2[0] - v1[0];
    ans[1] = v2[1] - v1[1];
    ans[2] = v2[2] - v1[2];
}

/// Determine the cross product
void VecMath::cross(real* ans, real* v1, real* v2) {
    ans[0] = v1[1]*v2[2] - v1[2]*v2[1];
    ans[1] = v1[2]*v2[0] - v1[0]*v2[2];
    ans[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

real VecMath::dot(real* v1, real* v2) {
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

/// Determine the l2 norm of a vector v
real VecMath::l2(real* v) {
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);    
}

/** Determine the area of a triangle using herons 
  * formula.
  * s = (a+b+c)/2
  * A = root(s(s-a)(s-b)(s-c))
  * \param p1 1st pt (real[3])
  * \param p2 2nd pt (real[3])
  * \param p3 3rd pt (real[3])
  * \returns The area  
  */
real VecMath::triArea(real* p1, real* p2, real* p3) {
    // This isn't particularly efficient, is it? 4 
    // sqrts?
    real a = dist(p1, p2);
    real b = dist(p2, p3);
    real c = dist(p3, p1);    
    real s = (a+b+c)*0.5;
    return sqrt(s*(s-a)*(s-b)*(s-c));   
}

