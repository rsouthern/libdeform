#ifndef FRAMEBUFFER_H
#define FRAMEBUFFER_H

#include <GL/glew.h>
#define GL_GLEXT_PROTOTYPES 
#include <GL/glext.h>
#include <GL/glut.h>

#include <string>

using namespace std;

/**
  * Encapsulation of the frame buffer
  */
class Framebuffer {
public:
    // Constructors
    Framebuffer();                  //< Construct an empty frame buffer
    ~Framebuffer();                 //< Destroy / delete the frame buffer

    // Main interface
    bool bind();                    //< Main interface - activate / draw to the buffer
    void unbind();                  //< Main interface - de-activate the buffer
    bool init(GLuint /*w*/, GLuint /*h*/); //< Construct a framebuffer with given width / height
    bool isInit() {return (fboId != 0);} //< Tell us if this is ready
    bool resize(GLuint, GLuint);    //< Resize the buffer

    // Check / print status
    bool checkFramebufferStatus();  //< Return true if we're good to go
    void printFramebufferInfo();    //< Dump the frame buffer stuff to the screen

    // The width and height of this particular frame buffer
    GLuint width() {return w;}
    GLuint height() {return h;}



private:
    GLuint fboId;                   //< The frame buffer id
    GLuint *rboId;                  //< The render buffers (color, depth and possibly stencil)

    GLuint w,h;                     //< The width and height of the FBO

    /// Internal functions for printing / checking the framebuffer
    string getRenderbufferParameters(GLuint id);
    string getTextureParameters(GLuint id);
    string convertInternalFormatToString(GLenum format);
};

#endif // FRAMEBUFFER_H


