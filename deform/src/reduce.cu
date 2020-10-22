#include "cuda_helper.cuh"
#include "reduce.cuh"
//#include <cutil_inline.h>
#include <cuda_runtime_api.h>
#include <helper_cuda.h>
#include <helper_functions.h>

extern uint NTHREADS;

/// An optimised version of the reduce (scan) function from GPU Gems 3.
template <unsigned int blockSize> __global__ void reduce(float *g_idata, float *g_odata, unsigned int n);

/**
*	Ripped from NVIDIA reduction example.
*
*  This version adds multiple elements per thread sequentially.  This reduces the overall
*  cost of the algorithm while keeping the work complexity O(n) and the step complexity O(log n).
* (Brent's Theorem optimization)
* \param g_idata (input) Vector to add up. Actual length must be 2^n.
* \param g_odata (input/output) The scratch space for adding using the parallel sum. The first entry contains the sum.
* \param n (input) The length of the vectors to add.
* \see reduceVector()
*/
template <unsigned int blockSize> __global__ void reduce(float *g_idata, float *g_odata, unsigned int n) {
    extern __shared__ float sdata[];

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(blockSize*2) + threadIdx.x;
    unsigned int gridSize = blockSize*2*gridDim.x;
    sdata[tid] = 0;

    // we reduce multiple elements per thread.  The number is determined by the 
    // number of active thread blocks (via gridSize).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
    while (i < n)
    {
        sdata[tid] += g_idata[i] + g_idata[i+blockSize];  
        i += gridSize;
    } 
    __syncthreads();

    // do reduction in shared mem
	if (blockSize >= 2048) { if (tid < 1024) { sdata[tid] += sdata[tid + 1024]; } __syncthreads(); }
	if (blockSize >= 1024) { if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
    if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
    if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
    if (blockSize >= 128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } __syncthreads(); }
    
#ifndef __DEVICE_EMULATION__
    if (tid < 32)
#endif
    {
        if (blockSize >=  64) { sdata[tid] += sdata[tid + 32]; EMUSYNC; }
        if (blockSize >=  32) { sdata[tid] += sdata[tid + 16]; EMUSYNC; }
        if (blockSize >=  16) { sdata[tid] += sdata[tid +  8]; EMUSYNC; }
        if (blockSize >=   8) { sdata[tid] += sdata[tid +  4]; EMUSYNC; }
        if (blockSize >=   4) { sdata[tid] += sdata[tid +  2]; EMUSYNC; }
        if (blockSize >=   2) { sdata[tid] += sdata[tid +  1]; EMUSYNC; }
    }
    
    // write result for this block to global mem 
    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

/**
  * Wrap the function which sums a vector together.
  * Calls a GPU accelerated sum function.
  * \param V Vector to sum.
  * \param N Length of the vector to sum
  * \returns The sum of vector V.
  */
extern "C" float reduceVector(float *V, const uint N) {
	// First deduce the optimum number of threads & blocks
	uint nThreads, nBlocks;
	
	if (N == 1) {
		// Determine the best number of threads here!
		nThreads = 1;
	} else {
		nThreads = (N < (NTHREADS*2)) ? (N/2) : NTHREADS;
	}
	nBlocks = N/(2*nThreads);

	int smemSize = nThreads * sizeof(float);
	float *work, *res, ret_val;

	// One piece of work for each block
	checkCudaErrors(cudaMalloc(&work, sizeof(float)*nBlocks));
	
	switch (nThreads) {
		case 2048:
			reduce<2048><<< nBlocks, nThreads, smemSize >>>(V, work, N); break;
		case 1024:
			reduce<1024><<< nBlocks, nThreads, smemSize >>>(V, work, N); break;
        case 512:
            reduce<512><<< nBlocks, nThreads, smemSize >>>(V, work, N); break;
        case 256:
            reduce<256><<< nBlocks, nThreads, smemSize >>>(V, work, N); break;
        case 128:
            reduce<128><<< nBlocks, nThreads, smemSize >>>(V, work, N); break;
        case 64:
            reduce< 64><<< nBlocks, nThreads, smemSize >>>(V, work, N); break;
        case 32:
            reduce< 32><<< nBlocks, nThreads, smemSize >>>(V, work, N); break;
        case 16:
            reduce< 16><<< nBlocks, nThreads, smemSize >>>(V, work, N); break;
        case  8:
            reduce<  8><<< nBlocks, nThreads, smemSize >>>(V, work, N); break;
        case  4:
            reduce<  4><<< nBlocks, nThreads, smemSize >>>(V, work, N); break;
        case  2:
            reduce<  2><<< nBlocks, nThreads, smemSize >>>(V, work, N); break;
        case  1:
            reduce<  1><<< nBlocks, nThreads, smemSize >>>(V, work, N); break;
        }
	res = (float*)malloc(sizeof(float)*nBlocks);
	checkCudaErrors(cudaMemcpy(res,work,sizeof(float)*nBlocks,cudaMemcpyDeviceToHost));
	ret_val = res[0];
	free(res);
	cudaFree(work);
	return ret_val;
}

