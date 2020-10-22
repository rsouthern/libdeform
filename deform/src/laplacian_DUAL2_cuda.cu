#include <laplacian_DUAL2_cuda.cuh>
#include <stdio.h>
//#include <cutil_math.cuh>
#include <cuda_helper.cuh>
#include <cuda_tools.cuh>
#include <cusp_tools.cuh>
#include <reduce.cuh>
#include <helper_cuda.h>
#include <helper_functions.h>
#include <helper_math.h>


//#include <cuPrintf.cu>

extern uint NTHREADS;

/// A device function to find the 4x4 determinant of matrix A.
inline __host__ __device__ float det4x4(const float A[]);

namespace CudaDual2 {

/// Update a single force vector
__global__ void correctOneFace_dual2(const uint nX,
                                     const uint nTri,
                                     const uint nB,
                                     const real_t *X,
                                     const uint *Tri,
                                     real_t *N0,
                                     real_t *fA,
                                     real_t *F_in,
                                     real_t *F_ext,
                                     real_t *B,
                                     real_t pressure = real_t(0.0));

/// Determine the volume of one face element
__global__ void faceVolume(const uint nX,
                           const uint nTri,
                           const uint *Tri,
                           const real_t *X,
                           float *Vol);


uint nX, nTri, nB;
real_t *X, *F_in, *F_ext, *B, *N0, *fA;
float *fVol; // Face volume
float V0;
uint *Tri;

/**
 * \param nX number of points, nrows of X, h and the modified rows of B
 * \param nTri number of faces, nrows of Tri, fA
 * \param Tri Triangle face data, size nTri*3
 * \param N0 Initial normal associated with each dual face
 * \param fA Area of each face, size nTri
 */
extern "C" void init_dual2_cu(const uint _nX,
                              const uint _nTri,
                              const uint _nB,
                              const uint *_Tri,
                              const double *_X,
                              const double *_F_in,
                              const double *_F_ext,
                              const double *_N0,
                              const double *_fA) {
    nX = _nX;
    nTri = _nTri;
    nB = _nB;

    // Allocate chunks of memory
    checkCudaErrors(cudaMalloc(&X, sizeof(real_t)*nX*3));
    checkCudaErrors(cudaMalloc(&F_in, sizeof(real_t)*nTri*3));
    checkCudaErrors(cudaMalloc(&F_ext, sizeof(real_t)*nTri*3));
    checkCudaErrors(cudaMalloc(&B, sizeof(real_t)*nB*3));
    checkCudaErrors(cudaMalloc(&Tri, sizeof(uint)*nTri*3));
    checkCudaErrors(cudaMalloc(&N0, sizeof(real_t)*nTri*3));
    checkCudaErrors(cudaMalloc(&fA, sizeof(real_t)*nTri));
    checkCudaErrors(cudaMalloc(&fVol, sizeof(float)*nTri));

    // Triangle data can be copied directly
    checkCudaErrors(cudaMemcpy(Tri, _Tri, sizeof(uint)*nTri*3,cudaMemcpyHostToDevice));

    // If we are not supporting double precision, the type must be cast (tedious)
    uint i,j;
    if (sizeof(real_t) != sizeof(double)) {
        real_t *f_N0 = (real_t*)malloc(sizeof(real_t)*nTri*3);
        real_t *f_fA = (real_t*)malloc(sizeof(real_t)*nTri);
        real_t *f_F_in = (real_t*)malloc(sizeof(real_t)*nTri*3);
        real_t *f_F_ext = (real_t*)malloc(sizeof(real_t)*nTri*3);
        real_t *f_X = (real_t*)malloc(sizeof(real_t)*nX*3);

        for (i=0; i<nTri; ++i) {
            for (j=0; j<3; ++j) {
                f_N0[i+nTri*j] = (real_t) _N0[i+nTri*j];
                f_F_in[i+nTri*j] = (real_t) _F_in[i+nTri*j];
                f_F_ext[i+nTri*j] = (real_t) _F_ext[i+nTri*j];
            }
            f_fA[i] = (real_t) _fA[i];
        }
        // Copy across the initial X as well
        for (i=0; i<nX; ++i)
            for (j=0; j<3; ++j)
                f_X[i+nX*j] = (real_t) _X[i+nX*j];
        checkCudaErrors(cudaMemcpy(N0, f_N0, sizeof(real_t)*nTri*3,cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(fA, f_fA, sizeof(real_t)*nTri, cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(F_in, f_F_in, sizeof(real_t)*nTri*3, cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(F_ext, f_F_ext, sizeof(real_t)*nTri*3, cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(X, f_X, sizeof(real_t)*nX*3, cudaMemcpyHostToDevice));
        free(f_N0);
        free(f_fA);
        free(f_F_in);
        free(f_F_ext);
        free(f_X);
    } else {
        checkCudaErrors(cudaMemcpy(N0, _N0, sizeof(real_t)*nTri*3,cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(fA, _fA, sizeof(real_t)*nTri, cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(F_in, _F_in, sizeof(real_t)*nTri*3, cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(F_ext, _F_ext, sizeof(real_t)*nTri*3, cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(X, _X, sizeof(real_t)*nX*3, cudaMemcpyHostToDevice));
    }
    // Compute the enclosed volume from the initial mesh state
    V0 = enclosedVolume();
    fprintf(stderr, "\nInitial Enclosed Volume = %f", V0);

    // Initialise the ability to printf from within the loop
    //checkCudaErrors(cudaPrintfInit());
}

/**
  * Clear away used memory
  */
extern "C" void shutdown_dual2_cu() {
    // Clear away memory
    checkCudaErrors(cudaFree(X));
    checkCudaErrors(cudaFree(F_in));
    checkCudaErrors(cudaFree(F_ext));
    checkCudaErrors(cudaFree(B));
    checkCudaErrors(cudaFree(Tri));
    checkCudaErrors(cudaFree(N0));
    checkCudaErrors(cudaFree(fA));
    checkCudaErrors(cudaFree(fVol));

    // Finish with CUDA device
    //cudaPrintfEnd();

    // Set it all to NULL
    X = F_in = F_ext = N0 = B = fA = NULL;
    Tri = NULL;
}

/**
 * Correct
 */
extern "C" void correct_dual2_cu(const double *_X, double *_B, double pressure_factor) {
    // Copy this all onto the GPU
    uint i,j;
    real_t *f_X, *f_B;
    if (sizeof(real_t) != sizeof(double)) {
        f_X = (real_t*) malloc(sizeof(real_t)*nX*3);
        f_B = (real_t*) malloc(sizeof(real_t)*nB*3);
        for (i=0; i<nX*3; ++i) f_X[i] = (real_t) _X[i];
        for (i=0; i<nB*3; ++j) f_B[i] = (real_t) _B[i];
        checkCudaErrors(cudaMemcpy(X, f_X, sizeof(real_t)*nX*3, cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(B, f_B, sizeof(real_t)*nB*3, cudaMemcpyHostToDevice));
        free(f_X); // Don't delete B - we need it later
    } else {
        checkCudaErrors(cudaMemcpy(X, _X, sizeof(real_t)*nX*3, cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(B, _B, sizeof(real_t)*nB*3, cudaMemcpyHostToDevice));
    }

    // Call the parallel function
    uint nThreads = 512;
    uint nBlocks = nTri / nThreads + 1;

    real_t A = totalArea();
    float V = enclosedVolume();
    real_t pressure = pressure_factor * real_t(V/V0 - 1.0f)/A;
    fprintf(stderr, "\nV0=%f, V=%f, A=%f, pressure=%f", V0, V, A, pressure);

    correctOneFace_dual2<<<nBlocks, nThreads>>>(nX, nTri, nB, X, Tri, N0, fA, F_in, F_ext, B, pressure);
    //checkCudaErrors(cudaPrintfDisplay(stderr));
    // Read B back from the GPU for use in the external application
    if (sizeof(real_t) != sizeof(double)) {
        checkCudaErrors(cudaMemcpy(f_B, B, sizeof(real_t)*nB*3, cudaMemcpyDeviceToHost));
        for (i=0; i<nX*3; ++i) {
            _B[i] = (double) f_B[i];
        }
        free(f_B);
    } else {
        checkCudaErrors(cudaMemcpy(_B, B, sizeof(real_t)*nB*3, cudaMemcpyDeviceToHost));
    }
    //cudaPrintfDisplay(stderr,true);
}

/**
  * Compute the total volume in the closed manifold
  */
extern "C" float enclosedVolume() {
    // Create a workspace in which to write our volume data
    float *h_fVol = (float*) malloc(sizeof(float)*nTri);

    // Compute the volume vector of local tetrahedral volumes (positive and negative!)
    uint nThreads = 512;
    uint nBlocks = nTri / nThreads + 1;
    faceVolume<<<nBlocks, nThreads>>>(nX, nTri, Tri, X, fVol);

    // Sum the results
    checkCudaErrors(cudaMemcpy(h_fVol, fVol, sizeof(float)*nTri, cudaMemcpyDeviceToHost));
    uint i;
    float ret_val = 0.0f;
    for (i=0; i<nTri; ++i) {
        ret_val += h_fVol[i];
    }
    //float ret_val = reduceVector(h_fVol, nTri);
    free(h_fVol);
    return ret_val;
}

/**
  * Determine total surface area
  */
extern "C" real_t totalArea() {
    real_t *h_fA = (real_t*) malloc(sizeof(real_t)*nTri);
    checkCudaErrors(cudaMemcpy(h_fA, fA, sizeof(real_t)*nTri, cudaMemcpyDeviceToHost));
    uint i;
    real_t ret_val = real_t(0.0);
    for (i=0; i<nTri; ++i) {
        ret_val += h_fA[i];
    }
    return ret_val;
}

/**
 * Correct
 */
__global__ void correctOneFace_dual2(const uint nX,
                                     const uint nTri,
                                     const uint nB,
                                     const real_t *X,
                                     const uint *Tri,
                                     real_t *N0,
                                     real_t *fA,
                                     real_t *F_in,
                                     real_t *F_ext,
                                     real_t *B,
                                     real_t pressure) {

    // Calculate the face index
    uint fidx = blockIdx.x*blockDim.x + threadIdx.x;

    // Check we haven't overrun our thread count
    if (fidx >= nTri) return;

    // Our old face normal - assume normalized!
    float3 old_normal = make_float3(N0[fidx+nTri*0], N0[fidx+nTri*1], N0[fidx+nTri*2]);

    // Compute the local face normal
    uint a = Tri[fidx+0*nTri];
    uint b = Tri[fidx+1*nTri];
    uint c = Tri[fidx+2*nTri];

    // Determine each of the points
    float3 v1 = make_float3(X[b+0*nX]-X[a+0*nX], X[b+1*nX]-X[a+1*nX] ,X[b+2*nX]-X[a+2*nX]);
    float3 v2 = make_float3(X[c+0*nX]-X[a+0*nX], X[c+1*nX]-X[a+1*nX], X[c+2*nX]-X[a+2*nX]);
    float3 normal = cross(v1,v2);

    //float3 normal = faceNormal(nX, X, nTri, Tri, fidx);
    //cuPrintf("\nNormal(%d)=[%f,%f,%f]", fidx, normal.x, normal.y, normal.z);


    // This normal has not been normalised - compute the magnitude and use this as the area
    float area = sqrtf(dot(normal,normal));

    // Slow normalisation - is it actually faster to use rsqrtf?
    normal = normal / area;

    // Determine the rotation angle
    float phi = acosf(dot(normal,old_normal));

    // Determine the axis of rotation and normalise
    float3 axis = cross(normal, old_normal);
    float inv_axis_norm = rsqrtf(axis.x*axis.x + axis.y*axis.y + axis.z*axis.z);
    axis = axis * inv_axis_norm;

    // Compute the rotation matrix
    float R[9];                        // A 3x3 matrix to store the rotation
    float3 f_in = make_float3(F_in[fidx+nTri*0], F_in[fidx+nTri*1], F_in[fidx+nTri*2]);
    float3 f_ext = make_float3(F_ext[fidx+nTri*0], F_ext[fidx+nTri*1], F_ext[fidx+nTri*2]);
    float3 res_in = f_in;
    float3 res_ext = f_ext;

    if (phi > 0.00001) {
        rotationMatrix(axis, phi, R);
        res_in = rotateVector(R, f_in);
    }

    // Determine the scaling factor    
    real_t scale = sqrtf(area/(fA[fidx]));
    //real_t scale = 1.0f;

    // Copy the result back to the output vector
    F_in[fidx+nTri*0] = scale * res_in.x;
    F_in[fidx+nTri*1] = scale * res_in.y;
    F_in[fidx+nTri*2] = scale * res_in.z;

    res_ext = (pressure*area) * normal ;
    F_ext[fidx+nTri*0] = res_ext.x;
    F_ext[fidx+nTri*1] = res_ext.y;
    F_ext[fidx+nTri*2] = res_ext.z;

    // Update the answer vector B
    uint i;
/*
    cuPrintf("\nB_before(%d)=[%f,%f,%f], B_after(%d)=[%f,%f,%f]",
             fidx, B[fidx+nB*0], B[fidx+nB*1], B[fidx+nB*2],
             fidx, F_in[fidx+nTri*0], F_in[fidx+nTri*1], F_in[fidx+nTri*2]);
*/
    for (i=0; i<3; ++i) {
        B[fidx+nB*i] = F_in[fidx+nTri*i] + F_ext[fidx+nTri*i];
    }

    // Update the normal direction
    N0[fidx+nTri*0] = normal.x;
    N0[fidx+nTri*1] = normal.y;
    N0[fidx+nTri*2] = normal.z;

    // Update the face area
    fA[fidx] = area;
}

/**
  * Compute the determinant of a 4x4 matrix A.
  */
inline __host__ __device__ float det4x4(const float A[]) {
    float k1 = A[10]*A[15] - A[11]*A[14];
    float k2 = A[6]*A[15] - A[7]*A[14];
    float k3 = A[6]*A[11] - A[7]*A[10];
    float k4 = A[2]*A[15] - A[14]*A[3];
    float k5 = A[2]*A[11] - A[10]*A[3];
    float k6 = A[2]*A[7] - A[3]*A[6];

    return A[0]*(A[5]*k1 - A[9]*k2 + A[13]*k3) -
           A[4]*(A[1]*k1 - A[9]*k4 + A[13]*k5) +
           A[8]*(A[1]*k2 - A[5]*k4 + A[13]*k6) -
           A[12]*(A[1]*k3 - A[5]*k5 + A[9]*k6);
}

/**
  * Compute the face volume of this particular face by creating a 4x4 matrix and finding it's determinant.
  */
__global__ void faceVolume(const uint nX, const uint nTri, const uint *Tri, const real_t *X, float *Vol) {
    uint fidx = blockIdx.x*blockDim.x + threadIdx.x;
    if (fidx >= nTri) return;

    // The default point from which to determine volume is set as [0,0,0].
    // The function det4x4 could be optimised accordingly.
    float A[] = {0.0f, 0.0f, 0.0f, 1.0f,
                 0.0f, 0.0f, 0.0f, 1.0f,
                 0.0f, 0.0f ,0.0f, 1.0f,
                 0.0f, 0.0f, 0.0f, 1.0f};


    uint i,j, vidx;
    for (i=0; i<3; ++i) {
        vidx = Tri[fidx+i*nTri];
        for (j=0; j<3; ++j) {
            if (vidx < nX) {
                A[j+i*4] = (float) X[vidx+j*nX];
            } else printf("ERROR!!!!!!! fidx=%d, vidx=%d", fidx, vidx);
        }
    }

    //printf("\nfaceVolume(%d): nX=%d, nTri=%d", fidx, nX, nTri);


    // Return the determinant of the 4X4 matrix as the local volume
    Vol[fidx] = det4x4(A);
    //printf("\n det([%f,%f,%f,%f;%f,%f,%f,%f;%f,%f,%f,%f;%f,%f,%f,%f]')=%f ",
    //       A[0],A[1],A[2],A[3],A[4],A[5],A[6],A[7],A[8],A[9],A[10],A[11],A[12],A[13],A[14],A[15],Vol[fidx]);
}

}

