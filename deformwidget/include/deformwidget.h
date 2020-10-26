#ifndef _DEFORMWIDGET_H
#define _DEFORMWIDGET_H

//== INCLUDES =================================================================

// GL includes
#include <GL/glew.h>
#include <GL/glext.h>
#include <GL/glut.h>

// QT includes
#include <QWidget>
#include <QString>
#include <QMessageBox>
#include <QFileDialog>
#include <QPushButton>
#include <QRadioButton>
#include <QButtonGroup>
#include <QVBoxLayout>

// OpenMesh includes
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Utils/getopt.h>
#include <OpenMesh/Tools/Utils/Timer.hh>
#include <OpenMesh/Apps/QtViewer/MeshViewerWidgetT.hh>

// Stdlib includes
#include <list>

// Include for cholmesh package
#include <mesh.h>

// Our own class for managing frame buffers
#include <framebuffer.h>
#include <vertexcolormap.h>

// A deformer class
#include <laplacian.h>

// A trackball class for deformation
#include <choltrackball.h>

// The names of the different Laplacian operators
#define DEF_LAPLACIAN       0
#define DEF_LAPLACIAN_COT   1
#define DEF_LAPLACIAN_DUAL  2
#define DEF_LAPLACIAN_DUAL2 3

using namespace std;
using namespace OpenMesh;  
using namespace OpenMesh::Attributes;

typedef MeshViewerWidgetT<Mesh> DeformWidgetParent;

// Define the key used for deformation
#define KEY_SELECT      Qt::Key_Shift
#define KEY_DEFORM      Qt::Key_Space

// Define the types of selections
#define STYPE_NOTHING     0
#define STYPE_ANCHORS     1
#define STYPE_CONSTRAINTS 2

// Application state
#define STATE_IDLE      0
#define STATE_SELECT    1
#define STATE_DEFORM    2

typedef struct SelectionBox {
    int x;
    int y;
    int radius;
    int type;
    SelectionBox() : x(0), y(0), radius(10), type(STYPE_NOTHING) {}
    void draw();
} SelectionBox;


/**
  *
  */
class DeformWidget : public DeformWidgetParent
{
    Q_OBJECT
public:
    /// default constructor
    DeformWidget(QWidget* parent=0);
    ~DeformWidget();

    OpenMesh::IO::Options& options() { return _options; }    
    const OpenMesh::IO::Options& options() const { return _options; }
    void setOptions(const OpenMesh::IO::Options& opts) { _options = opts; }

    void open_mesh_gui(QString fname);
    void open_texture_gui(QString fname);

    /// Override the open mesh function
    virtual bool open_mesh(const char* /*fname*/, OpenMesh::IO::Options /*_opt*/);

    /// Override the OpenMesh render function for our special drawing routine. Called when updateGL is called in the parent.
    virtual void draw_scene(const std::string & /*_drawmode*/);

    /// Called whenever resizeGL is called in the parent
    virtual void resize_scene(const int &_w, const int &_h);

    virtual void draw_openmesh(const std::string& /*_drawmode*/);

public slots:
    void query_open_mesh_file();
    void query_open_texture_file();
    void set_deformer();
    void clear_deformer();
    void set_front_face(bool);

protected:
    // A custom event filter
    virtual bool eventFilter(QObject *, QEvent *);

    vector<bool> constraints;       //< Constraints are user handles
    vector<bool> anchors;           //< Anchors don't move (much)

private:
    OpenMesh::IO::Options _options;
    uint state;                     //< Either IDLE, SELECT or DEFORM
    uint currentDeformer;           //< Keep track of the current deformer mode
    bool selectionFrontFace;        //< True if you can only select front facing vertices
    Framebuffer fbo;                //< Contains a frame buffer object
    SelectionBox sbox;              //< A selfcontained selection thingy
    bool fboInit;                   //< Keep track of whether we've initialised the frame buffer
    unsigned char *unique_vertex_colors;   //< Array storing our unique vertex colors
    VertexColorMap cmap;            //< A color map from a vertex to an index and visa versa
    Laplacian *op;                  //< Our local Laplacian deformer operator
    void renderSelectionImage();    //< Render a special image for the purposes of vertex / face selection
    void generate_unique_vertex_colors(); //< Generate unique vertex colors for selection
    void find_selected_vertices(uint, uint, uint, vector<bool> *); //< Find vertices in the selection radius


    QPointF pixelPosToViewPos(const QPointF&); //< Transform the position for the trackball
    CholTrackBall deformTrackBall;  //< A trackball to handle deformation of anchors
};


#endif //_DEFORM_WIDGET
