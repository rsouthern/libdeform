#ifndef __DEFINES_H
#define __DEFINES_H

/**
  * Unsigned int type - conflicts with others?
  */
typedef unsigned int uint;

/**
  * Precision to use for calculations. Note that this must be changed to float when
  * CUDA compute capability < 1.3. This should be part of cuda_helper and cuda functions
  * should use templates. One day I will do this. Note the type real conflicts with others.
  */
typedef double real_t;


/**
  * Output some more info from our function
  */
#define VERBOSE

/**
  * The maximum number of threads to use for sanity sake
  */
#define MAX_THREADS 1024

/**
  * Debugging command if you want to save the big matrices to an ASCII format for debugging
  */
#define SAVE_MATRICES

/**
  * The null hash indicates the point isn't in the grid (this shouldn't happen if your extents are correctly chosen)
  */
#define NULL_HASH			UINT_MAX

/**
  * Syncthreads doesn't work in emulation mode, so disable it if necessary
  */
#ifdef __DEVICE_EMULATION__
#define EMUSYNC __syncthreads()
#else
#define EMUSYNC
#endif

/**
 * The default tolerance to determine convergence.
 */
#define DEFAULT_TOL 1.0e-4f

/**
  * The default number of iterations.
  */ 
#define DEFAULT_MAX_ITER 1000

#ifndef NULL
#define NULL 0
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#endif

