#ifndef LAPLACIAN_DUAL2_CUDA_CUH
#define LAPLACIAN_DUAL2_CUDA_CUH

//#include <cutil_inline.h>
#include <defines.h>

namespace CudaDual2 {

/// Initialise the memory needed for the problem
extern "C" void init_dual2_cu(const uint _nX,
                              const uint _nTri,
                              const uint _nB,
                              const uint *_Tri,                              
                              const double *_X,
                              const double *_F_in,
                              const double *_F_ext,
                              const double *_N0,
                              const double *_fA);

/// Shutdown cuda and delete the allocated memory
extern "C" void shutdown_dual2_cu();

/// Correct the force vectors in F_in and F_ext
extern "C" void correct_dual2_cu(const double */*_X*/, double */*B*/, double /*pressure_factor*/ = 0.0);

/// Compute the current volume enclosed by the mesh
extern "C" float enclosedVolume();

/// Sum each face area to get the total
extern "C" real_t totalArea();


}




#endif // LAPLACIAN_DUAL2_CUDA_CUH


