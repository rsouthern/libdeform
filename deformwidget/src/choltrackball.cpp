#include "choltrackball.h"
#include <stdio.h>

CholTrackBall::CholTrackBall() : glTrackball() {
    init();
}

/**
  * Initialise the trackball
  */
void CholTrackBall::init() {
    // Get the common cholmod workspace
    c = CholmodSharedSingleton::Instance()->getCholmodCommon();

    // Create a rotation matrix and set initial values
    Mat = cholmod_allocate_dense(4,4,4,CHOLMOD_REAL,c);

    // A good idea to do this!
    trans[0]=trans[1]=trans[2]=0.0f;
}

/**
  * Delete allocated memory and remove reference to cholmod
  */
CholTrackBall::~CholTrackBall() {
    cholmod_free_dense(&Mat, c);
    CholmodSharedSingleton::Instance()->stop();
}

void CholTrackBall::setGLView() {
    // Call the parent function first
    glTrackball::setGLView();

    // Now apply the panning transform
    glTranslated(trans[0],trans[1],trans[2]);
}


/**
  * Our only contribution: Return the rotation as a cholmod matrix
  */
cholmod_dense *CholTrackBall::transformMat() {
    double x = result.v[0];
    double y = result.v[1];
    double z = result.v[2];
    double w = result.s;

    double *Mx = (double*) Mat->x;
    Mx[0] = 1.0-2.0*y*y-2.0*z*z;
    Mx[4] = 2.0*x*y-2.0*z*w;
    Mx[8] = 2.0*x*z + 2.0*y*w;
    Mx[12] = trans[0];

    Mx[1] = 2.0*x*y + 2.0*z*w;
    Mx[5] = 1.0-2.0*x*x-2.0*z*z;
    Mx[9] = 2.0*y*z - 2.0*z*w;
    Mx[13] = trans[1];

    Mx[2] = 2.0*x*z - 2.0*y*w;
    Mx[6] = 2.0*y*z + 2.0*x*w;
    Mx[10] = 1.0 - 2.0*x*x - 2.0*y*y;
    Mx[14] = trans[2];

    Mx[3] = Mx[7] = Mx[11] = 0.0;
    Mx[15] = 1.0;
    return Mat;
}


void CholTrackBall::MouseMove(int x, int y) {
    glTrackball::MouseMove(x,y);
    lastx = x;
    lasty = y;
    return;

    switch (dragMode) {
    case dmPan:
        // This requires the modelview and projection to be on the stack already

        // Retrieve our matrices
        GLdouble model[16], proj[16];
        GLint view[4];
        glGetDoublev(GL_MODELVIEW_MATRIX, model);
        glGetDoublev(GL_PROJECTION_MATRIX, proj);
        glGetIntegerv(GL_VIEWPORT, view);

        // Now unproject our last point and this point
        GLdouble startPos[3];
        GLdouble endPos[3];
        GLdouble diff[3];
        GLdouble winZ;
        int i;

        winZ = 0.99;
        //glReadPixels(lastx, view[3]-lasty, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);
        //winZ = 1.0 - (0.01 / trackBall.getModelDistance());
        gluUnProject((GLdouble) lastx, view[3]-lasty, winZ,
                     model, proj, view,
                     &(startPos[0]),&(startPos[1]),&(startPos[2]));

        //glReadPixels(x, view[3]-y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);
        gluUnProject((GLdouble) x, view[3]-y, winZ,
                     model, proj, view,
                     &(endPos[0]),&(endPos[1]),&(endPos[2]));
        for (i = 0; i < 3; ++i)
            trans[i] = endPos[i] - startPos[i];
        break;
    default:
        glTrackball::MouseMove(x,y);
        return;
    }

    /* update the last x and y */
    lastx = x;
    lasty = y;
}

void CholTrackBall::reset() {
    trans[0]=trans[1]=trans[2]=0.0;
    glTrackball::reset();
}

