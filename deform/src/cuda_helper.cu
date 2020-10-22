//#include <cutil_math.cuh>
#include <cuda_runtime_api.h>
#include <cuda_helper.cuh>
#include <cuda_tools.cuh>
#include <helper_cuda.h>
#include <helper_functions.h>

// Define a global variable stashing the number of threads to use
uint NTHREADS;

#ifdef USE_CULA
/**
  * Useful function to check the results of a call to CULA.
  */
extern "C" void culaCheckStatus(culaStatus status) {
	if(!status)
		return;
	printf("CULA Failure: ");
	if(status == culaArgumentError)
		printf("Invalid value for parameter %d\n", culaGetErrorInfo());
	else if(status == culaDataError)
		printf("Data error (%d)\n", culaGetErrorInfo());
	else if(status == culaBlasError)
		printf("Blas error (%d)\n", culaGetErrorInfo());
	else if(status == culaRuntimeError)
		printf("Runtime error (%d)\n", culaGetErrorInfo());
	else
		printf("%s\n", culaGetStatusString(status));

	//getchar();
	culaShutdown();
	exit(EXIT_FAILURE);
}
#endif //USE_CULA

/**
  * Initialise the CUDA device. Currently this initialises the best device in the system, but doesn't worry about
  * shared memory. One day this may be an issue!
  * \return The best number of threads to use.
  */
extern "C" int initCudaDevice() {
	// Query the hardware
	cudaDeviceProp prop;
	// Use the first device by default
	int nDevices = 0;
    checkCudaErrors(cudaGetDeviceCount(&nDevices));
	int i;
	int bestDevice = -1;
	int bestThreads = 0;
	for (i=0; i < nDevices; ++i) {
        checkCudaErrors(cudaGetDeviceProperties(&prop, 0));
		printf("\nDevice %d: %s, regsPerBlock=%d, sharedMemPerBlock=%d, maxThreadsPerBlock=%d",
                        (int) i, prop.name, (int) prop.regsPerBlock, (int) prop.sharedMemPerBlock,  (int) prop.maxThreadsPerBlock);
		if ((bestDevice == -1) || (prop.maxThreadsPerBlock > bestThreads)) {
			bestDevice = i;
			bestThreads = prop.maxThreadsPerBlock;
		}
	}
	if (bestDevice == -1) {
		printf("\ninitCudaDevice(): Fatal error - no CUDA device found!\n");
		return -1;
	} else {
		printf("\ninitCudaDevice(): Using device %d.\n", bestDevice);
	}
    checkCudaErrors(cudaSetDevice(bestDevice));
	if (bestThreads > MAX_THREADS) {
		NTHREADS = MAX_THREADS;
	} else {
		NTHREADS = bestThreads;
	}

#ifdef USE_CULA
	culaStatus status;
	status = culaInitialize();
	if(status != culaNoError) {
		printf("\ninitCudaDevice(): Error initialising CULA. Try undefining USE_CULA.\n");
		return status;
	}
#endif //USE_CULA
	return 0;
}

/**
  * Shutdown the CUDA device and disable CULA if necessary.
  */
extern "C" int shutdownCudaDevice() {
#ifdef USE_CULA
	culaShutdown();
#endif //USE_CULA	
	return 0;
}

