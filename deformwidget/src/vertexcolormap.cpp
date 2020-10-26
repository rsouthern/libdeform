#include "vertexcolormap.h"
#include <stdio.h>
#include <cmath>
#include <math.h>
#include <limits>

void CheckOpenGLError(const char* stmt, const char* fname, int line)
{
    GLenum err = glGetError();
    if (err != GL_NO_ERROR)
    {
        printf("OpenGL error %08x, at %s:%i - for %s\n", err, fname, line, stmt);
        abort();
    }
}

#ifdef _DEBUG
    #define GL_CHECK(stmt) do { \
            stmt; \
            CheckOpenGLError(#stmt, __FILE__, __LINE__); \
        } while (0)
#else
    #define GL_CHECK(stmt) stmt
#endif

VertexColorMap::VertexColorMap() {
    // Retrieve the color depth
    GLint redBits, greenBits, blueBits;

    // These calls appear to be deprecated and return unreliable data. Setting this to 8-bit as default.
    GL_CHECK(glGetIntegerv (GL_RED_BITS, &redBits));
    GL_CHECK(glGetIntegerv (GL_GREEN_BITS, &greenBits));
    GL_CHECK(glGetIntegerv (GL_BLUE_BITS, &blueBits));
    redBits = greenBits = blueBits = 8;

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

string VertexColorMap::toBinary(const GLuint &input) {
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
bool VertexColorMap::colorToIdx(const unsigned char &r, const unsigned char &g, const unsigned char &b, GLuint *idx) {
    if (r == 0 && g == 0 && b == 0) return false;
    (*idx) = (r  << redShift) + (g  << greenShift) + (b << blueShift) - 1;
    return true;
}

// A map from the index to a color
bool VertexColorMap::idxToColor(const GLuint &_idx, unsigned char *r, unsigned char *g, unsigned char *b) {
    // Always add one to the input and remove one from the output
    GLuint idx = 1 + _idx;
    (*r) = ((idx) & redMask) >> redShift;
    (*g) = ((idx) & greenMask) >> greenShift;
    (*b) = ((idx) & blueMask) >> blueShift;
    return true;
}

GLuint VertexColorMap::makeMask(const GLuint &_x) {
    GLuint res = GLuint(0);
    GLuint x = _x;
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

