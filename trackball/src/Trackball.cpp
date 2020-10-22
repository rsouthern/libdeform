#include "Trackball.h"

Trackball::Trackball()
{
    trackballSize = 0.8;
    result = gmQuaternion(1.0,0.0,0.0,0.0);

}

gmQuaternion Trackball::move(float p1x, float p1y, float p2x, float p2y)
{
    gmVector3 a;            /* Axis of rotation */
    gmVector3 p1, p2, d;
    gmQuaternion temp;

    float phi;              /* how much to rotate about axis */
    float t;

    if (p1x == p2x && p1y == p2y)
       return gmQuaternion(1.0,0.0,0.0,0.0);

    /* First, figure out z-coordinates for projection of P1 
	  * and P2 to deformed sphere.  
	  */
    p1 = gmVector3(p1x,p1y,projectToSphere(p1x,p1y));
    p2 = gmVector3(p2x,p2y,projectToSphere(p2x,p2y));

    a = cross(p2,p1);
    d = p1 - p2;
    t = d.length()/(2.0*trackballSize);
    t = CLAMP(t,-1.0,1.0); 
    phi = 2.0 * asin(t);
    temp.FromAngleAxis(phi,a);
    result = result * temp;
    return result;
}

void Trackball::ToRotationMatrix (gmMatrix4 &R) const
{
  result.ToRotationMatrix (R);
}

float Trackball::projectToSphere(float x, float y)
{
    float d, t, z;

    d = sqrt(x*x + y*y);
    if (d < trackballSize * 0.70710678118654752440) {    /* Inside sphere */
        z = sqrt(trackballSize*trackballSize - d*d);
    } else {                                 /* On hyperbola */
        t = trackballSize / 1.41421356237309504880;
        z = t*t / d;
    }
    return z;
}

void Trackball::setTrackballSize(float radius)
{
   trackballSize = radius;
}

void Trackball::reset() {
	result = gmQuaternion(1.0,0.0,0.0,0.0);
}
