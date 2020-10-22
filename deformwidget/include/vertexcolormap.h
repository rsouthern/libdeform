#ifndef VERTEXCOLORMAP_H
#define VERTEXCOLORMAP_H

#include <GL/glew.h>
#include <GL/glext.h>
#include <GL/glut.h>

#include <stdio.h>
#include <string>

using namespace std;

/**
  *
  */
class VertexColorMap {
public:
    VertexColorMap();
    ~VertexColorMap();
    bool colorToIdx(unsigned char /*r*/, unsigned char /*g*/, unsigned char /*b*/, GLuint */*idx*/);
    bool idxToColor(GLuint /*idx*/, unsigned char */*r*/, unsigned char */*g*/, unsigned char */*b*/);
private:
    GLuint redMask, greenMask, blueMask;
    GLuint makeMask(GLuint);
    string toBinary(GLuint /*input*/);      //< For debugging
    int redShift, greenShift, blueShift;
};


#endif // VERTEXCOLORMAP_H
