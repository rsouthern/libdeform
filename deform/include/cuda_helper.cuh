#ifndef _CUDA_HELPER_H
#define _CUDA_HELPER_H

#include <defines.h>

#ifdef USE_CULA
#include <cula.h>
extern "C" void culaCheckStatus(culaStatus status);
#endif //USE_CULA

/// Initialise your best cuda device and return the maximum number of supported threads
extern "C" int initCudaDevice();

/// Shutdown the device - currently only matters if you are using CULA
extern "C" int shutdownCudaDevice();

#endif //_CUDA_HELPER_H
