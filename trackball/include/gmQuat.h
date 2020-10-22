// Extended libgm++ routine - Graphics Math Library
// Modified from Magic-Sofware.com code.
// Caleb Lyness
// August 1999

#ifndef GMQUAT_H_CAL
#define GMQUAT_H_CAL
#undef NDEBUG

#include "gmUtils.h"
#include "gmVec3.h"
#include "gmMat3.h"
#include "gmMat4.h"

class gmQuaternion {
  //protected:
  public:
    gmVector3 v;	/* vector part of quaternion */
    double    s;	/* scalar part */

  public:
    gmQuaternion ();
    gmQuaternion (double s, double x, double y, double z);
    gmQuaternion (double s, const gmVector3 &v);
    gmQuaternion (const gmQuaternion& q);

    /* conversion between quaternions, matrices, and angle-axes */
    void FromRotationMatrix (const gmMatrix4& R);
    void ToRotationMatrix (gmMatrix4 &R) const;
    void FromAngleAxis (const double& angle, const double& ax, const double& ay, const double& az);
    void FromAngleAxis (const double angle, const gmVector3& point);
    void ToAngleAxis (double& angle, double& ax, double& ay, double& az) const;
    void ToAngleAxis (double& angle, gmVector3& point) const;

    /* arithmetic operations */
    gmQuaternion operator+ (const gmQuaternion& q) const;
    gmQuaternion operator- (const gmQuaternion& q) const;
    gmQuaternion operator* (const gmQuaternion& q) const;
    gmQuaternion operator/ (gmQuaternion q) const; /* multiply by the inverse! */
    gmQuaternion operator* (double c) const;
    gmQuaternion operator/ (double c) const;
    friend gmQuaternion operator* (double c, const gmQuaternion& q);
    gmQuaternion operator- () const;

    /* assignment and aritmetic operations */
    gmQuaternion& operator= (const gmQuaternion& q);
    gmQuaternion& operator+= (const gmQuaternion& q);
    gmQuaternion& operator-= (const gmQuaternion& q);
    gmQuaternion& operator*= (const gmQuaternion& q);
    gmQuaternion& operator/= (const gmQuaternion& q);
    gmQuaternion& operator*= (double c);
    gmQuaternion& operator/= (double c);

    /* functions of a quaternion */
    double Dot (const gmQuaternion& q) const;  // dot product of quaternions
    double Norm () const;  // squared-length of quaternion
    gmQuaternion Inverse () const;  // apply to non-zero quaternion
    gmQuaternion UnitInverse () const;  // apply to unit-length quaternion
    gmQuaternion Exp () const;
    gmQuaternion Log () const;

	 void normalize();

    /* rotation of a point by a quaternion */
    gmVector3 operator* (const gmVector3& pt) const;

    /* spherical linear interpolation */
    static gmQuaternion Slerp (double t, const gmQuaternion& p,
        const gmQuaternion& q);

    static gmQuaternion SlerpExtraSpins (double t, const gmQuaternion& p,
        const gmQuaternion& q, int iExtraSpins);

    /* setup for spherical quadratic interpolation */
    static void Intermediate (const gmQuaternion& q0, const gmQuaternion& q1,
        const gmQuaternion& q2, gmQuaternion& a, gmQuaternion& b);

    /* spherical quadratic interpolation */
    static gmQuaternion Squad (double t, const gmQuaternion& p,
        const gmQuaternion& a, const gmQuaternion& b, const gmQuaternion& q);

    /* output */
    friend std::ostream & operator << ( std::ostream &, const gmQuaternion&);
};

class gmQuatInterp
{
  public:
    gmQuatInterp();
};

// ------- Define the inline routines: ---------------------------
inline 
std::ostream& operator << (std::ostream& os, const gmQuaternion& q)
{
  os << "[ " << q.s << "," << q.v << " ]";
  return os;
}

inline 
gmQuaternion::gmQuaternion()
  :v(0.0,0.0,0.0) { s = 0.0; }

inline
gmQuaternion::gmQuaternion (double s, const gmVector3 &v)
{
  gmQuaternion::v = v;
  gmQuaternion::s = s;
}

inline
gmQuaternion::gmQuaternion (double s, double x, double y, double z)
  :v(x,y,z)
{
  gmQuaternion::s = s;
}

inline
gmQuaternion::gmQuaternion (const gmQuaternion& q)
{
  v = q.v;
  s = q.s;
}

inline void 
gmQuaternion::FromAngleAxis (const double angle, const gmVector3& point)
{ FromAngleAxis(angle,point[0],point[1],point[2]); }

inline void 
gmQuaternion::ToAngleAxis (double& angle, gmVector3& point) const
{ ToAngleAxis(angle,point[0],point[1],point[2]); }

inline gmQuaternion
gmQuaternion::operator+ (const gmQuaternion& q) const
{ return (gmQuaternion(s+q.s,v+q.v)); }

inline gmQuaternion
gmQuaternion::operator- (const gmQuaternion& q) const
{ return (gmQuaternion(s-q.s,v-q.v)); }

inline gmQuaternion
gmQuaternion::operator* (double c) const
{ return (gmQuaternion(c*s,c*v)); }

inline gmQuaternion
gmQuaternion::operator/ (gmQuaternion q) const /* multiply by the inverse! */
{ return ((*this)*q.Inverse()); }

inline gmQuaternion
gmQuaternion::operator/ (double c) const
{ return ((*this)*(1/c)); }

inline
gmQuaternion operator* (double c, const gmQuaternion& q)
{ return (q*c); }

inline gmQuaternion
gmQuaternion::operator- () const
{ return gmQuaternion(-s,-v); }

inline gmQuaternion&
gmQuaternion::operator= (const gmQuaternion& q)
{ v=q.v; s=q.s; return *this; }

inline gmQuaternion&
gmQuaternion::operator+= (const gmQuaternion& q)
{ v+=q.v; s+=q.s; return *this; }

inline gmQuaternion&
gmQuaternion::operator-= (const gmQuaternion& q)
{ v-=q.v; s-=q.s; return *this; }

inline gmQuaternion&
gmQuaternion::operator*= (const gmQuaternion& q)
{ return (*this = (*this * q));}

inline gmQuaternion&
gmQuaternion::operator/= (const gmQuaternion& q)
{ return (*this = *this/q);}

inline gmQuaternion&
gmQuaternion::operator*= (double c)
{ return (*this = *this*c); }

inline gmQuaternion&
gmQuaternion::operator/= (double c)
{ return (*this = *this * (1/c));}

/* functions of a quaternion */
inline double 
gmQuaternion::Dot (const gmQuaternion& q) const
{ return (s*q.s + dot(v,q.v)); }

inline double 
gmQuaternion::Norm () const
{ return (gmSqr(s)+v.lengthSquared()); }

inline gmQuaternion 
gmQuaternion::UnitInverse () const  // assert:  'this' is unit length
{  return gmQuaternion(s,-v); }

#endif
