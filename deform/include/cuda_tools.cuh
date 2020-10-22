#ifndef CUDA_TOOLS_CUH
#define CUDA_TOOLS_CUH

//#include <cutil_inline.h>
//#include <cutil_math.h>
#include <defines.h>
#include <helper_math.h>

/*
inline __host__ __device__ float3 cross(float3 v1, float3 v2) {
    return make_float3(v1.y*v2.z-v1.z*v2.y,
                       v1.z*v2.x-v1.x*v2.z,
                       v1.x*v2.y-v1.y*v2.x);
}*/

/**
  * Compute the face normal
  */
inline __host__ __device__ float3 faceNormal(const uint nX, const real_t *X, const uint nTri, const uint *Tri, const int _idx) {
    float3 normal = make_float3(0.0f,0.0f,0.0f);
    uint idx;
    if (_idx > nTri) return normal;

    // If the index is the default value, deduce the index from the current thread
    if (_idx < 0) idx = blockIdx.x*blockDim.x + threadIdx.x;

    uint a = Tri[idx+0*nTri];
    uint b = Tri[idx+1*nTri];
    uint c = Tri[idx+2*nTri];

    // Determine each of the points    
    float3 v1 = make_float3(X[b+0*nX]-X[a+0*nX], X[b+1*nX]-X[a+1*nX] ,X[b+2*nX]-X[a+2*nX]);
    float3 v2 = make_float3(X[c+0*nX]-X[a+0*nX], X[c+1*nX]-X[a+1*nX], X[c+2*nX]-X[a+2*nX]);
    normal = cross(v1,v2);
    return normal;
}

/**
  * Construct the rotation matrix out of an up vector and an angle of
  * rotation using eulers method
  * http://mathworld.wolfram.com/EulerParameters.html.
  */
inline __host__ __device__ void rotationMatrix(const float3 N, const float phi, float *R) {
    float sin_phi = sin(0.5*phi);
    float e[] = {cos(0.5*phi), N.x*sin_phi, N.y*sin_phi, N.z*sin_phi};
    float ee[] = {e[0]*e[0], e[1]*e[1], e[2]*e[2], e[3]*e[3]};

    R[0] = ee[0] + ee[1] - ee[2] - ee[3];
    R[1] = real_t(2)*(e[1]*e[2] - e[0]*e[3]);
    R[2] = real_t(2)*(e[1]*e[3] + e[0]*e[2]);
    R[3] = real_t(2)*(e[1]*e[2] + e[0]*e[3]);
    R[4] = ee[0] - ee[1] + ee[2] - ee[3];
    R[5] = real_t(2)*(e[2]*e[3] - e[0]*e[1]);
    R[6] = real_t(2)*(e[1]*e[3] - e[0]*e[2]);
    R[7] = real_t(2)*(e[2]*e[3] + e[0]*e[1]);
    R[8] = ee[0] - ee[1] - ee[2] + ee[3];
}

/**
  * Rotate a vector by the given rotation matrix
  */
inline __host__ __device__ float3 rotateVector(const float *R, const float3 vec) {    
    return make_float3(R[0]*vec.x + R[3]*vec.y + R[6]*vec.z,
                       R[1]*vec.x + R[4]*vec.y + R[7]*vec.z,
                       R[2]*vec.x + R[5]*vec.y + R[8]*vec.z);
}


#endif //CUDA_TOOLS_CUH


