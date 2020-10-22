#include "laplacian_DUAL_cuda.cuh"
#include <stdio.h>
#include <cuda_helper.cuh>
#include <cuda_tools.cuh>
//#include "/usr/local/cuda/sdk/C/src/simplePrintf/cuPrintf.cu"
#include "cusp_tools.cuh"
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>
#include <helper_cuda_gl.h>
#include <helper_math.h>

namespace CudaDual {

/// Correct a single vector thread in B
__global__ void correctOneFace_dual(const uint nX, const uint nTri, const uint nB, const real_t *X,
                                    const uint *Tri, const real_t *h, real_t *fA, real_t *B);

uint nX, nTri, nB;
real_t *X, *B, *h, *fA;
uint *Tri;

/*
extern "C" void test_cusp(cholmod_factor *fac) {
    dMat dM;
    if (cholmod_factor_to_cusp(fac, &dM)) {
        printf("\nSuccess!");
    }
}*/

/*
 * \param nX number of points, nrows of X, h and the modified rows of B
 * \param nTri number of faces, nrows of Tri, fA
 * \param nB number of rows in B, including constraints (only the first nX will be modified)
 * \param X Point data stored in column major order, size nX*3
 * \param Tri Triangle face data, size nTri*3
 * \param h height of each feature, size nTri
 * \param fA Area of each face, size nTri
 * \param B Detail vector data which will be updated as a result of this function
 */

extern "C" void init_dual_cu(const uint _nX,
                                const uint _nB,
                                const uint _nTri,
                                const uint *_Tri,
                                const double *_h,
                                const double *_fA) {
    nX = _nX;
    nTri = _nTri;
    nB = _nB;
    checkCudaErrors(cudaMalloc(&X, sizeof(real_t)*nX*3));
    checkCudaErrors(cudaMalloc(&B, sizeof(real_t)*nB*3));
    checkCudaErrors(cudaMalloc(&Tri, sizeof(uint)*nTri*3));
    checkCudaErrors(cudaMalloc(&h, sizeof(real_t)*nTri));
    checkCudaErrors(cudaMalloc(&fA, sizeof(real_t)*nTri));

    // Triangle data can be copied directly
    checkCudaErrors(cudaMemcpy(Tri, _Tri, sizeof(uint)*nTri*3,cudaMemcpyHostToDevice));

    // If we are not supporting double precision, the type must be cast (tedious)
    uint i;
    if (sizeof(real_t) != sizeof(double)) {
        real_t *f_h = (real_t*)malloc(sizeof(real_t)*nTri);
        real_t *f_fA = (real_t*)malloc(sizeof(real_t)*nTri);
        for (i=0; i<nTri; ++i) {
            f_h[i] = (real_t) _h[i];
            f_fA[i] = (real_t) _fA[i];
        }
        checkCudaErrors(cudaMemcpy(h, f_h, sizeof(real_t)*nTri,cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(fA, f_fA, sizeof(real_t)*nTri, cudaMemcpyHostToDevice));
        free(f_h);
        free(f_fA);
    } else {
        checkCudaErrors(cudaMemcpy(h, _h, sizeof(real_t)*nTri, cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(fA, _fA, sizeof(real_t)*nTri, cudaMemcpyHostToDevice));
    }

    //cudaPrintfInit();
}

/**
  * Clear away used memory
  */
extern "C" void shutdown_dual_cu() {
    // Clear away memory
    checkCudaErrors(cudaFree(X));
    checkCudaErrors(cudaFree(B));
    checkCudaErrors(cudaFree(Tri));
    checkCudaErrors(cudaFree(h));
    checkCudaErrors(cudaFree(fA));

    // Finish with CUDA device
    //cudaPrintfEnd();

    // Set it all to NULL
    X = B = h = fA = NULL;
    Tri = NULL;
}

/**
 * Correct the matrix B by assigning to it the value -h_i * n_i.
 */
extern "C" void correct_dual_cu(const double *_X, double *_B) {
    //fprintf(stderr, "\ncorrect_cuda()");
    // Copy this all onto the GPU

    uint i;
    real_t *f_X, *f_B;
    if (sizeof(real_t) != sizeof(double)) {
        f_X = (real_t*)malloc(sizeof(real_t)*nX*3);
        f_B = (real_t*)malloc(sizeof(real_t)*nB*3);
        for (i=0; i<nX*3; ++i) f_X[i] = (real_t) _X[i];
        for (i=0; i<nB*3; ++i) f_B[i] = (real_t) _B[i];
        checkCudaErrors(cudaMemcpy(X, f_X, sizeof(real_t)*nX*3, cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(B, f_B, sizeof(real_t)*nB*3, cudaMemcpyHostToDevice));
        free(f_X);
    } else {
        checkCudaErrors(cudaMemcpy(X, _X, sizeof(real_t)*nX*3, cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(B, _B, sizeof(real_t)*nB*3, cudaMemcpyHostToDevice));
    }

    // Call the parallel function
    uint nThreads = 512;
    uint nBlocks = nTri / nThreads + 1;
    correctOneFace_dual<<<nBlocks, nThreads>>>(nX, nTri, nB, X, Tri, h, fA, B);

    // Read B back from the GPU for use in the external application
    if (sizeof(real_t) != sizeof(double)) {
        checkCudaErrors(cudaMemcpy(f_B, B, sizeof(real_t)*nX*3, cudaMemcpyDeviceToHost));
        for (i=0; i<nB*3; ++i) _B[i] = (double) f_B[i];
        free(f_B);
    } else {
       checkCudaErrors(cudaMemcpy(_B, B, sizeof(float)*nB*3, cudaMemcpyDeviceToHost));
    }
}

/**
  * For each face in the mesh, compute the face normal and recompute the local feature vector.
  */
__global__ void correctOneFace_dual(const uint nX,
                             const uint nTri,
                             const uint nB,
                             const real_t *X,
                             const uint *Tri,
                             const real_t *h,
                             real_t *fA,
                             real_t *B) {
    // Calculate the face index
    uint fidx = blockIdx.x*blockDim.x + threadIdx.x;

    // Check we haven't overrun our thread count
    if (fidx >= nTri) return;

    // Compute the local face normal
    uint a = Tri[fidx+0*nTri];
    uint b = Tri[fidx+1*nTri];
    uint c = Tri[fidx+2*nTri];

    // Determine each of the points
    float3 v1 = make_float3(X[b+0*nX]-X[a+0*nX], X[b+1*nX]-X[a+1*nX] ,X[b+2*nX]-X[a+2*nX]);
    float3 v2 = make_float3(X[c+0*nX]-X[a+0*nX], X[c+1*nX]-X[a+1*nX], X[c+2*nX]-X[a+2*nX]);
    float3 normal = cross(v1,v2);

    // This normal has not been normalised - compute the magnitude and use this as the area
    float area = sqrtf(dot(normal,normal));

    // Slow normalisation - is it actually faster to use rsqrtf?
    normal = normal / area;

    // Determine the scaling factor
    real_t scale = sqrtf(area/(fA[fidx]));
    //real_t scale = 1.0f;

    // Scale our output vectors
    B[fidx+0*nB] = -scale * h[fidx] * ((real_t)normal.x);
    B[fidx+1*nB] = -scale * h[fidx] * ((real_t)normal.y);
    B[fidx+2*nB] = -scale * h[fidx] * ((real_t)normal.z);

    // Update the face area
    fA[fidx] = area;

    //cuPrintf("\ncorrectFaces(%d) - normal [%f,%f,%f], scale=%f, height=%f", fidx, normal.x, normal.y, normal.z, scale, h[fidx]);
}

}
