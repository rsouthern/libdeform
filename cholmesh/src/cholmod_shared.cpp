#include "cholmod_shared.h"


/// Global static pointer used to ensure a single instance of the class.
CholmodSharedSingleton* CholmodSharedSingleton::_Instance = NULL;

/// Initialise our reference counter to zero
unsigned int CholmodSharedSingleton::refCnt = 0;

/** This function is called to create an instance of the class.
    Calling the constructor publicly is not allowed. The constructor
    is private and is only called by this Instance function.
*/
CholmodSharedSingleton* CholmodSharedSingleton::Instance() {
    // Only allow one instance of class to be generated.
    if (!_Instance) _Instance = new CholmodSharedSingleton;
    ++refCnt;
    return _Instance;
}

/**
  * Stop using this singleton class. Needs to be explicitely called by the client.
  */
void CholmodSharedSingleton::stop() {
    // Sanity check - this could be a major problem!
    if (refCnt == 0) {
        fprintf(stderr, "\nCholmodSharedSingleton::stop() - Warning: stop called too many times on this singleton!");
        return;
    }
    --refCnt;

    if (refCnt == 0) {
        delete _Instance;
        _Instance = NULL;
    }
}

/**
  * Private constructor which will instantiate a cholmod_common structure and start cholmod.
  *
  * Note: Only DOUBLE precision is currently supported by CHOLMOD! This means that there
  * will be additional overhead typecasting floats to doubles when transfering CHOLMOD data
  * to and from the GPU. On Fermi cores this may no longer be a problem.
  */
CholmodSharedSingleton::CholmodSharedSingleton() {
    cholmod_start(&c);
    c.metis_memory = 4.0;
    c.itype = CHOLMOD_INT;
    c.dtype = CHOLMOD_DOUBLE; // Double precision! (make sure this corresponds to the definition of real in defines.h)
    c.print = 4;
}

/**
  * Private destructor to shutdown cholmod
  */
CholmodSharedSingleton::~CholmodSharedSingleton() {
    cholmod_finish(&c);
}


