#ifndef _LAPLACIAN_DUAL_CUDA_CUH
#define _LAPLACIAN_DUAL_CUDA_CUH

//#include <cutil_inline.h>
#include <defines.h>
#include <cholmod_shared.h>
#include <reduce.cuh>

namespace CudaDual {

/// Allocate memory
extern "C" void init_dual_cu(const uint _nX,
                                const uint _nB,
                                const uint _nTri,
                                const uint *_Tri,
                                const real_t *_h,
                                const real_t *_fA);

/// Shutdown by deleting memory
extern "C" void shutdown_dual_cu();

/// Given the current point data X and the answer vector B, update B
extern "C" void correct_dual_cu(const double *X, double *_B);

}

/// A testing function
extern "C" void test_cusp(cholmod_factor *fac);

#endif //_LAPLACIAN_DUAL_CUDA_CUH
