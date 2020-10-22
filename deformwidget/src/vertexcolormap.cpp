#include "vertexcolormap.h"
#include <stdio.h>
#include <cmath>
#include <math.h>
#include <limits>

VertexColorMap::VertexColorMap() {
    // Retrieve the color depth
    GLint redBits, greenBits, blueBits;
    glGetIntegerv (GL_RED_BITS, &redBits);
    glGetIntegerv (GL_GREEN_BITS, &greenBits);
    glGetIntegerv (GL_BLUE_BITS, &blueBits);

    fprintf(stderr, "\nRGB bits: %d,%d,%d", redBits,greenBits,blueBits);

    redMask = makeMask(redBits) << (greenBits + blueBits);
    greenMask = makeMask(greenBits) << blueBits;
    blueMask = makeMask(blueBits);

    // create a string from each mask
    redShift = (greenBits+blueBits);
    greenShift = (blueBits);
    blueShift = 0;

    fprintf(stderr, "\nVertexColorMap() - rgbBits=[%d,%d,%d], rgbMask=[%s,%s,%s], rgbShift=[%d,%d,%d]",
            redBits, greenBits, blueBits,
            toBinary(redMask).c_str(),
            toBinary(greenMask).c_str(),
            toBinary(blueMask).c_str(),
            redShift,greenShift,blueShift);
}

VertexColorMap::~VertexColorMap() {
}

string VertexColorMap::toBinary(GLuint input) {
    if (input == 0) return "0"; // trivial case
    string result;
    int i;
    for (i = numeric_limits<GLuint>::digits - 1; i >= 0; --i) {
        if (input & (1 << i)) result += "1";
        else if(!result.empty()) result += "0";
    }
    return result;
}

// A map from the color to the index
bool VertexColorMap::colorToIdx(unsigned char r, unsigned char g, unsigned char b, GLuint *idx) {
    if (r == 0 && g == 0 && b == 0) return false;
    (*idx) = (r  << redShift) + (g  << greenShift) + (b << blueShift) - 1;
    return true;
}

// A map from the index to a color
bool VertexColorMap::idxToColor(GLuint _idx, unsigned char *r, unsigned char *g, unsigned char *b) {
    // Always add one to the input and remove one from the output
    GLuint idx = 1 + _idx;
    (*r) = ((idx) & redMask) >> redShift;
    (*g) = ((idx) & greenMask) >> greenShift;
    (*b) = ((idx) & blueMask) >> blueShift;
    return true;
}

GLuint VertexColorMap::makeMask(GLuint x) {
    GLuint res = GLuint(0);
    if (x == 0) {
        fprintf(stderr, "makemask() called with zero value \n");
        return 0;
    }
    while (x > 0) {
        x--;
        res = (res << 1) | 0x01;
    }
    return res;
}

