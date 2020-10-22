#ifndef _trackball_h_CAL_
#define _trackball_h_CAL_

/*
 * Trackball code:
 *
 * Implementation of a virtual trackball.
 * Implemented by Gavin Bell, lots of ideas from Thant Tessman and
 *   the August '88 issue of Siggraph's "Computer Graphics," pp. 121-129.
 */

/* Ammended and updated for ease of use by Dennis Burford 07/11/00
*/

#include "gm.h"
#include "gmQuat.h"
#include "common.h"

class Trackball 
{
  public:
     Trackball();
	 virtual ~Trackball() {};
	 
	  gmQuaternion move(float p1x, float p1y, float p2x, float p2y);
    void ToRotationMatrix (gmMatrix4 &R) const;

	  void setTrackballSize(float radius);

	  virtual void reset();

  protected:
     gmQuaternion result;

  private:
     float trackballSize;

     float projectToSphere(float x, float y);
};

#endif
