#ifndef CHOLMOD_SHARED_H
#define CHOLMOD_SHARED_H

#include <cholmod.h>
#include <string>

/**
  * Usage:
  * In your constructor: cholmod_common *c = CholmodSharedSingleton::Instance()->getCholmodCommon;
  * In your destructor: CholmodSharedSingleton::stop();
  */
class CholmodSharedSingleton {
public:
   /// Call this to get our instance of the class. Constructs if necessary.
   static CholmodSharedSingleton* Instance();

   /// Call this to retrieve our single cholmod workspace
   cholmod_common* getCholmodCommon() {return &c;}

   /// Stop using this class (call this once in your destructor or something
   static void stop();

private:
   /// Cholmod common object
   cholmod_common c;

   /// Our one instance of the singleton
   static CholmodSharedSingleton* _Instance;

   /// Our reference counter
   static unsigned int refCnt;

   // Private members for our singleton
   CholmodSharedSingleton();                                              //< Private constructor
   ~CholmodSharedSingleton();                                             //< Private destructor(?)
   CholmodSharedSingleton(const CholmodSharedSingleton&);                 //< Prevent copy-construction
   CholmodSharedSingleton& operator=(const CholmodSharedSingleton&);      //< Prevent assignment
};

#endif // CHOLMOD_SHARED_H


