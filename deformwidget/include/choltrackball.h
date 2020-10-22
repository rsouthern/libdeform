#ifndef CHOLTRACKBALL_H
#define CHOLTRACKBALL_H

#include <cholmod_shared.h>
#include <glTrackball.h>


class CholTrackBall : public glTrackball {
public:
    CholTrackBall();
    ~CholTrackBall();

    /// Our only contribution: Return the rotation as a cholmod matrix
    cholmod_dense *rotationMat() const;

    /// Set the current GL view
    virtual void setGLView();

    /// Implement panning
    virtual void MouseMove(int x, int y);

    /// Retrieve a cholmod matrix for this transform
    cholmod_dense *transformMat();

    /// Override this one
    virtual void reset();

private:
    cholmod_dense *Mat;         //< Local storage for Cholmod transformation matrix
    cholmod_common *c;          //< Local pointer to the common structure
    GLfloat trans[3];           //< The translation stored for panning
    void init();                //< Initialise the cholmod matrix
};

#endif // CHOLTRACKBALL_H
