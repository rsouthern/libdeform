//---------------------------------------------------------------------------
#ifndef glTrackballH
#define glTrackballH

#include <GL/glut.h>
#include "Trackball.h"

enum tDragMode {
  dmTrackball,    /* Trackball mode: Left button down    */
  dmZoom,         /* Trackball mode: Middle button down  */
  dmPan,          /* RS: Pan mode: Right button down     */
  dmTranslate,    /* Fly Mode: Middle button down        */
  dmRotate,       /* Fly Mode: Left button down          */
  dmNone          /* No button down      */
};

//---------------------------------------------------------------------------
class glTrackball : public Trackball
{
  public:

    glTrackball();
	virtual ~glTrackball() {};
	
    virtual void setGLView();

    void setModelOrigin(float x, float y, float z);
    void setModelOrigin(gmVector3 o);

    void setModelDistance(float d);

    void setWindowDimensions(int w, int h);

    virtual void MouseMove(int x, int y);
    virtual void MouseButton(int button, int state, int x, int y);
    virtual void reset() {Trackball::reset();};
    float getModelDistance() {return modelDistance;};
    void draw();
    
  protected:
    int renorm_count;

    tDragMode dragMode;
    int lastx, lasty;

    gmVector3 modelOrigin;
    float modelDistance;

    int winw, winh;
};

//---------------------------------------------------------------------------
#endif
