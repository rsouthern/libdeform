//---------------------------------------------------------------------------
#include "glTrackball.h"
//---------------------------------------------------------------------------

glTrackball::glTrackball() : Trackball() {
  dragMode = dmNone;
  renorm_count = 97;
  modelOrigin = gmVector3(0.5, 0.5, 0.5);
  modelDistance = 1.0;
  winw = 512;
  winh = 512;
}


/**
  * Draw the trackball into the current render context
  */
void glTrackball::draw() {
    fprintf(stderr, "\nglTrackball::draw() - modelDistance=%f, modelOrigin=[%f,%f,%f]", modelDistance, modelOrigin[0],modelOrigin[1],modelOrigin[2]);

    glMatrixMode(GL_MODELVIEW);    
    glPushMatrix();

    setGLView();
    glColor3f(1.0f,1.0f,1.0f);
    glutWireSphere(modelDistance, 10,10);

    glColor3f(1.0f, 0.0f, 0.0f);
    glBegin(GL_LINE);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(modelDistance, 0.0f, 0.0f);
    glEnd(); //GL_LINE

    glColor3f(0.0f, 1.0f, 0.0f);
    glBegin(GL_LINE);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, modelDistance, 0.0f);
    glEnd(); //GL_LINE

    glColor3f(0.0f, 0.0f, 1.0f);
    glBegin(GL_LINE);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, modelDistance);
    glEnd(); //GL_LINE

    glPopMatrix();

}

//---------------------------------------------------------------------------

void glTrackball::setGLView()
{
   gmMatrix4 rotationMatrix;

   /* setup the view matrix */
   glTranslatef(0,0,-modelDistance);
   glTranslated(modelOrigin[0],modelOrigin[1],modelOrigin[2]);
   ToRotationMatrix(rotationMatrix);
   glMultMatrixd(rotationMatrix[0]);

   /* shift the model to the center of the trackball */
   glTranslated(-modelOrigin[0],-modelOrigin[1],-modelOrigin[2]);
}

//---------------------------------------------------------------------------

void glTrackball::setModelOrigin(float x, float y, float z)
{
  modelOrigin = gmVector3(x, y, z);
}

//---------------------------------------------------------------------------

void glTrackball::setModelOrigin(gmVector3 o)
{
  modelOrigin = o;
}

//---------------------------------------------------------------------------

void glTrackball::setModelDistance(float d)
{
  modelDistance = d;
}

//---------------------------------------------------------------------------

void glTrackball::setWindowDimensions(int w, int h)
{
  winw = w;
  winh = h;
}

//---------------------------------------------------------------------------

void glTrackball::MouseMove(int x, int y) {
   int len;
   static int count=0;
   GLfloat lastZoom;
   GLfloat x1,x2,y1,y2;

   switch (dragMode)
   {
      case dmTrackball:
      	   x1 = (GLfloat)(2*lastx-winw)/(GLfloat)winw;
           y1 = (GLfloat)(winh-2*lasty)/(GLfloat)winh;
           x2 = (GLfloat)(2*x-winw)/(GLfloat)winw;
           y2 = (GLfloat)(winh-2*y)/(GLfloat)winh;

           move(x1,y1,x2,y2);

           if (++count > renorm_count)
              /* Renomalise to prevent accum errors */
              result.normalize();
           break;

      case dmZoom:  /* "Z" motion */
           len = lasty - y;
           modelDistance += len*0.05;

		   /*
			// FOV zooming
           lastZoom = zoom;
           zoom += ((float)len/(float)h)*40.0;
           zoom = CLAMP(zoom,5.0,120.0);
           redisplay = (lastZoom != zoom);
		   //cout << quat << endl;
		  */
           break;
      default: return;
   }

   /* update the last x and y */
   lastx = x;
   lasty = y;
}

//---------------------------------------------------------------------------

void glTrackball::MouseButton(int button, int state, int x, int y)
{
   switch (button) {
     case 0: dragMode = (state == 1)? dmTrackball: dmNone; break; // left
     case 1: dragMode = (state == 1)? dmZoom: dmNone; break;      // middle
     case 2: dragMode = (state == 1)? dmPan: dmNone; break;       // right
     default: return;
   }

   /* update the last x and y */
   lastx = x;
   lasty = y;
}

//---------------------------------------------------------------------------

//#pragma package(smart_init)

