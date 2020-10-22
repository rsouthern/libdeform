#ifdef _MSC_VER
#  pragma warning(disable: 4267 4311)
#endif

#ifdef ARCH_DARWIN
#include <glew.h>
#include <glut.h>
#else
#define GL_GLEXT_PROTOTYPES
#include <GL/glew.h>
#include <GL/glext.h>
#include <GL/glut.h>
#endif

//#include <iostream>
//#include <fstream>
#include <QApplication>
#include <QMessageBox>
#include <QMainWindow>
#include <QMenuBar>
#include <QFileDialog>


#include "deformwidget.h"

#include "cuda_helper.cuh"
  
void create_menu(QMainWindow &w);
void usage_and_exit(int xcode);

int main(int argc, char **argv)
{
  // OpenGL check
//  QApplication::setColorSpec( QApplication::CustomColor );
  QApplication app(argc,argv);
#if !defined(__APPLE__)
  glutInit(&argc,argv);
  glewInit();
#endif


  if ( !QGLFormat::hasOpenGL() ) {
    QString msg = "System has no OpenGL support!";
    QMessageBox::critical( 0, QString("OpenGL"), msg + QString(argv[1]) );
    return -1;
  }

  int c;
  OpenMesh::IO::Options opt;
  
  while ( (c=getopt(argc,argv,"hbs"))!=-1 )
  {
     switch(c)
     {
       case 'b': opt += OpenMesh::IO::Options::Binary; break;
       case 'h':
          usage_and_exit(0);
       case 's': opt += OpenMesh::IO::Options::Swap; break;
       default:
          usage_and_exit(1);
     }
  }
  // create widget
  QMainWindow mainWin;
  DeformWidget w(&mainWin);
  w.setOptions(opt);
  mainWin.setCentralWidget(&w);

  create_menu(mainWin);

  // static mesh, hence use strips
  w.enable_strips();

  mainWin.resize(640, 480);
  mainWin.show(); 

  // load scene if specified on the command line
  if ( optind < argc )  
  {
    w.open_mesh_gui(argv[optind]);
  }

  if ( ++optind < argc )
  {
      w.open_texture_gui(argv[optind]);
  }

  // Start cuda (hate doing it here, but I get runtime errors when I try using a singleton class - apparently
  // memory is not being deallocated cleanly so we cannot call initCudaDevice() more than once!)
  if (initCudaDevice() != 0) {
      fprintf(stderr,"\nFailure to initialise CUDA!");
  }

  return app.exec();
}

void create_menu(QMainWindow &w)
{
    // Create the File menu
    using namespace Qt;
    QMenu *fileMenu = w.menuBar()->addMenu(w.tr("&File"));

    QAction* openAct = new QAction(w.tr("&Open mesh..."), &w);
    openAct->setShortcut(w.tr("Ctrl+O"));
    openAct->setStatusTip(w.tr("Open a mesh file"));
    QObject::connect(openAct, SIGNAL(triggered()), w.centralWidget(), SLOT(query_open_mesh_file()));
    fileMenu->addAction(openAct);

    QAction* texAct = new QAction(w.tr("Open &texture..."), &w);
    texAct->setShortcut(w.tr("Ctrl+T"));
    texAct->setStatusTip(w.tr("Open a texture file"));
    QObject::connect(texAct, SIGNAL(triggered()), w.centralWidget(), SLOT(query_open_texture_file()));
    fileMenu->addAction(texAct);

    QMenu *editMenu = w.menuBar()->addMenu(w.tr("&Edit"));
    QAction *setDeformerAct = new QAction(w.tr("&Set deformer"),&w);
    QAction *clearDeformerAct = new QAction(w.tr("&Clear anchors and constraints"),&w);
    QAction *setFrontFaceAct = new QAction(w.tr("&Front face selection"),&w);
    setFrontFaceAct->setCheckable(true);
    setFrontFaceAct->setChecked(true);
    QObject::connect(setDeformerAct, SIGNAL(triggered()), w.centralWidget(), SLOT(set_deformer()));
    QObject::connect(clearDeformerAct, SIGNAL(triggered()), w.centralWidget(), SLOT(clear_deformer()));
    QObject::connect(setFrontFaceAct, SIGNAL(toggled(bool)), w.centralWidget(), SLOT(set_front_face(bool)));

    editMenu->addAction(setDeformerAct);
    editMenu->addAction(clearDeformerAct);
    editMenu->addAction(setFrontFaceAct);
}

void usage_and_exit(int xcode)
{
   std::cout << "Usage: meshviewer [-s] [mesh] [texture]\n" << std::endl;
   std::cout << "Options:\n"
	     << "  -b\n"
	     << "    Assume input to be binary.\n\n"
             << "  -s\n"
             << "    Reverse byte order, when reading binary files.\n"
             << std::endl;
   exit(xcode);
}
