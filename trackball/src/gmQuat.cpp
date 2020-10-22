// gmQuaternion.cpp - quaternion class
//
// Extended libgm++ routine - Graphics Math Library
// Modified from Magic-Sofware.com code.
// Caleb Lyness
// August 1999

#include <math.h>
#include "gmQuat.h"

static double g_dEpsilon = 1e-03;  // cutoff for sin(angle) near zero
static double g_dPi = 4.0*atan(1.0);

void gmQuaternion::FromRotationMatrix (const gmMatrix4& R)
{
    double trace, root;
    int i,j,k;
    static int next[3] = { 1, 2, 0 };
    trace = R[0][0]+R[1][1]+R[2][2];

    if ( trace > 0.0 )
    {
        // |w| > 1/2, may as well choose w > 1/2
        root = sqrt(trace+1.0);  // 2w
        s = 0.5*root;
        root = 0.5/root;  // 1/(4w)
        v[0] = (R[2][1]-R[1][2])*root;
        v[1] = (R[0][2]-R[2][0])*root;
        v[2] = (R[1][0]-R[0][1])*root;
    }
    else
    {
        // |w| <= 1/2
	i = (R[1][1] > R[0][0])? 1: 0;
        if ( R[2][2] > R[i][i] ) i=2;
        j = next[i];
        k = next[j];

        root = sqrt(R[i][i]-R[j][j]-R[k][k]+1.0);
	v[i] = 0.5*root;
        root = 0.5/root;
        s = (R[k][j]-R[j][k])*root;
        v[j] = (R[j][i]+R[i][j])*root;
        v[k] = (R[k][i]+R[i][k])*root;
    }
}

void gmQuaternion::ToRotationMatrix (gmMatrix4 &R) const
{
    double tx  = 2.0*v[0];
    double ty  = 2.0*v[1];
    double tz  = 2.0*v[2];
    double twx = tx*s;
    double twy = ty*s;
    double twz = tz*s;
    double txx = tx*v[0];
    double txy = ty*v[0];
    double txz = tz*v[0];
    double tyy = ty*v[1];
    double tyz = tz*v[1];
    double tzz = tz*v[2];

    R[0][0] = 1.0-(tyy+tzz);
    R[0][1] = txy-twz;
    R[0][2] = txz+twy;
    R[0][3] = 0.0;

    R[1][0] = txy+twz;
    R[1][1] = 1.0-(txx+tzz);
    R[1][2] = tyz-twx;
    R[1][3] = 0.0;

    R[2][0] = txz-twy;
    R[2][1] = tyz+twx;
    R[2][2] = 1.0-(txx+tyy);
    R[2][3] = 0.0;

    R[3][0] = 0.0;
    R[3][1] = 0.0;
    R[3][2] = 0.0;
    R[3][3] = 1.0;
}

void gmQuaternion::FromAngleAxis (const double& angle, const double& ax,
                                  const double& ay, const double& az)
{
    // The quaternion representing the rotation is
    //  q = cos(A/2)+sin(A/2)*(x*i+y*j+z*k)

    double halfAngle = 0.5*angle;
    double sn = sin(halfAngle);
    double norm = 1/sqrt(ay*ay + ax*ax + az*az);

    v = gmVector3(sn*ax*norm,sn*ay*norm,sn*az*norm);
    s = cos(halfAngle);
}

void gmQuaternion::ToAngleAxis (double& angle, double& ax, double& ay,
                              double& az) const
{
        // The quaternion representing the rotation is
        //   q = cos(A/2)+sin(A/2)*(x*i+y*j+z*k)

    double length2 = v.lengthSquared();
    if ( length2 > 0.0 )
    {
        angle = 2.0*acos(s);
        double invlen = 1.0/sqrt(length2);
        ax = v[0]*invlen;
        ay = v[1]*invlen;
        az = v[2]*invlen;
    }
    else
    {
        // angle is 0 (mod 2*pi), so any axis will do
        angle = 0;
        ax = 1.0;
        ay = 0.0;
        az = 0.0;
    }
}

/*inline*/ gmQuaternion
gmQuaternion::operator* (const gmQuaternion& q) const
{
    
    return gmQuaternion
    (
       s*q.s - v[0]*q.v[0] - v[1]*q.v[1] - v[2]*q.v[2],
       s*q.v[0] + v[0]*q.s + v[1]*q.v[2] - v[2]*q.v[1],
       s*q.v[1] + v[1]*q.s + v[2]*q.v[0] - v[0]*q.v[2],
       s*q.v[2] + v[2]*q.s + v[0]*q.v[1] - v[1]*q.v[0]
    );
}

inline gmQuaternion 
gmQuaternion::Inverse () const
{
    double norm = Norm();
    assert (norm > 0.0);
    norm = 1.0/norm;
    return gmQuaternion(s*norm,-norm*v);
}

inline gmQuaternion 
gmQuaternion::Exp () const
{
    // If q = A*(x*i+y*j+z*k) where (x,y,z) is unit length, then
    // exp(q) = cos(A)+sin(A)*(x*i+y*j+z*k).  If sin(A) is near zero,
    // use exp(q) = cos(A)+A*(x*i+y*j+z*k) since A/sin(A) has limit 1.

    double angle = v.length();;
    double sn = sin(angle);
    return gmQuaternion(cos(angle), 
           ((fabs(sn)>=g_dEpsilon)? v*(sn/angle): v));
}

inline gmQuaternion
gmQuaternion::Log () const
{
    // If q = cos(A)+sin(A)*(x*i+y*j+z*k) where (x,y,z) is unit length, then
    // log(q) = A*(x*i+y*j+z*k).  If sin(A) is near zero, use log(q) =
    // sin(A)*(x*i+y*j+z*k) since sin(A)/A has limit 1.

    gmQuaternion result; /* s is by default 0.0 */

    if ( fabs(s) < 1.0 )
    {
        double angle = acos(s);
        double sn = sin(angle);
        if ( fabs(sn) >= g_dEpsilon )
        {
	    result.v = v*(angle/sn);
            return result;
        }
    }

    result.v = v;
    return result;
}

void 
gmQuaternion::normalize()
{
	float val;
	val = Norm();
   s /= val;	
	v /= val;
}

gmVector3 
gmQuaternion::operator* (const gmVector3& pt) const
{
    // Given a vector u = (x0,y0,z0) and a unit length quaternion
    // q = <w,x,y,z>, the vector v = (x1,y1,z1) which represents the
    // rotation of u by q is v = q*u*q^{-1} where * indicates quaternion
    // multiplication and where u is treated as the quaternion <0,x0,y0,z0>.
    // Note that q^{-1} = <w,-x,-y,-z>, so no real work is required to
    // invert q.  Now
    //
    //   q*u*q^{-1} = q*<0,x0,y0,z0>*q^{-1}
    //     = q*(x0*i+y0*j+z0*k)*q^{-1}
    //     = x0*(q*i*q^{-1})+y0*(q*j*q^{-1})+z0*(q*k*q^{-1})
    //
    // As 3-vectors, q*i*q^{-1}, q*j*q^{-1}, and 2*k*q^{-1} are the columns
    // of the rotation matrix computed in Quaternion::ToRotationMatrix.  The
    // vector v is obtained as the product of that rotation matrix with
    // vector u.  As such, the quaternion representation of a rotation
    // matrix requires less space than the matrix and more time to compute
    // the rotated vector.  Typical space-time tradeoff...

    gmMatrix4 R;
    ToRotationMatrix(R);

    return gmVector3(
        R[0][0]*pt[0] + R[0][1]*pt[1] + R[0][2]*pt[2], // x
        R[1][0]*pt[0] + R[1][1]*pt[1] + R[1][2]*pt[2], // y
        R[2][0]*pt[0] + R[2][1]*pt[1] + R[2][2]*pt[2]  // z
    );
}

gmQuaternion 
gmQuaternion::Slerp (double t, const gmQuaternion& p, const gmQuaternion& q)
{
    // assert:  p.Dot(q) >= 0 (guaranteed in NiRotKey::Interpolate methods)
    double dCos = p.Dot(q);

    // numerical round-off error could create problems in call to acos
    if ( dCos < -1.0 )
        dCos = -1.0;
    else if ( dCos > 1.0 )
        dCos = 1.0;

    double dAngle = acos(dCos);
    double dSin = sin(dAngle);  // fSin >= 0 since fCos >= 0

    if ( dSin < g_dEpsilon )
    {
        return p;
    }
    else
    {
        double dInvSin = 1.0/dSin;
        double dCoeff0 = sin((1.0-t)*dAngle)*dInvSin;
        double dCoeff1 = sin(t*dAngle)*dInvSin;
        return dCoeff0*p + dCoeff1*q;
    }
}

gmQuaternion 
gmQuaternion::SlerpExtraSpins (double t, const gmQuaternion& p,
    const gmQuaternion& q, int iExtraSpins)
{
    // assert:  p.Dot(q) >= 0 (guaranteed in RotationKey::Preprocess)
    double dCos = p.Dot(q);

    // numerical round-off error could create problems in call to acos
    if ( dCos < -1.0 )
        dCos = -1.0;
    else if ( dCos > 1.0 )
        dCos = 1.0;

    double dAngle = acos(dCos);
    double dSin = sin(dAngle);  // fSin >= 0 since fCos >= 0

    if ( dSin < g_dEpsilon )
    {
        return p;
    }
    else
    {
        double dPhase = g_dPi*iExtraSpins*t;
        double dInvSin = 1.0/dSin;
        double dCoeff0 = sin((1.0-t)*dAngle - dPhase)*dInvSin;
        double dCoeff1 = sin(t*dAngle + dPhase)*dInvSin;
        return dCoeff0*p + dCoeff1*q;
    }
}

void 
gmQuaternion::Intermediate (const gmQuaternion& q0, const gmQuaternion& q1,
    const gmQuaternion& q2, gmQuaternion& a, gmQuaternion& b)
{
    // assert:  q0, q1, q2 are unit quaternions

    gmQuaternion q0inv = q0.UnitInverse();
    gmQuaternion q1inv = q1.UnitInverse();
    gmQuaternion p0 = q0inv*q1;
    gmQuaternion p1 = q1inv*q2;
    gmQuaternion arg = 0.25*(p0.Log()-p1.Log());
    gmQuaternion marg = -arg;

    a = q1*arg.Exp();
    b = q1*marg.Exp();
}

gmQuaternion 
gmQuaternion::Squad (double t, const gmQuaternion& p,
    const gmQuaternion& a, const gmQuaternion& b, const gmQuaternion& q)
{
    return Slerp(2*t*(1-t),Slerp(t,p,q),Slerp(t,a,b));
}

