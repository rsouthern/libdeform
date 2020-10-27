# libdeform
### Laplacian mesh deformation library

This application provides 3 libraries and a Qt widget (based on the OpenMesh viewer) to visualise the libraries at work.
- *libcholmesh*: A link between OpenMesh (a half-edge mesh data structure) and Cholmod (a Cholesky factorisation package). Also supports dual mesh representations. Clever and potentially reusable.
- *libdeform*: Implementations of base Laplacian mesh deformation with either even or Cotangent weights, Dual Graph Laplacian (https://core.ac.uk/download/pdf/193853802.pdf) and 2nd Order Dual Graph Laplacian (as yet unpublished research code of limited utility). Solves the LU problem on the CPU, but updates vertex positions in realtime using CUDA. Object Orientation means that this *should* be portable to ARAP, ACAP or other modern day laplacian representations. Clever reusable code. 
- *libtrackball*: An ancient trackball library I inherited from Denis Burford in the UCT days. Uses deprecated OpenGL and desparately needs updating. Avoid reusing this please!
- *deformwidget*: A widget for visualising the final application based on the OpenMesh Viewer Widget (included in the packages dir and slightly modified). Uses deprecated OpenGL and is full of bugs. Desperately needs a modern port. Depends on the trackball heavily. Avoid reusing this please!

Note that there are decent implementations in libIGL of these operations (https://libigl.github.io/) but these are not GPU accelerated.

Compilation instructions:
- git clone <library path on git>
- cd libdeform
- mkdir build
- cd build
- cmake ..
- (assuming this works) make

Libs will be installed in <approot>/libs and the (only) binary will be in your <approot>/bin directory.
You will need an obj mesh to load.

Dependencies:
- OpenMesh (https://www.graphics.rwth-aachen.de/software/openmesh/download/) - define a OPENMESH_PATH environment variable pointing to the base directory.
- svdlibc (https://github.com/lucasmaystre/svdlibc) - define a SVDLIB_PATH environment variable pointing the base directory.
- SuiteSparse (https://github.com/DrTimothyAldenDavis/SuiteSparse/releases) - define a SUITESPARSE_PATH environment variable pointing to the base directory.
- CUDA (this was originally developed with CUDA 2.0)
- cuda-samples (https://github.com/NVIDIA/cuda-samples) - define a CUDA_SAMPLES_PATH environment variable pointing to the base directory.
- OpenGL, GLEW, GLU

Usage instructions:
- run <approot>/bin/deformwidget
- load a mesh from the menu
- select deformation type from the menu
- hold shift to select anchors and handles:
  - handles (in red) move when you deform and are constrained to follow the path of deformation. Hold middle button to select these.
  - anchors (in blue) "hold" the mesh in place. Hold left mouse button to select these. You should make a ring of these around your handles to localise deformation.
  - resize the selection box by scrolling the mouse wheel.
  - note that selecting nothing will reset handles (there is also a menu option to do this)
- hold space to perform the deformation:
  - click and drag the left mouse while space is depressed (oh I'm so depressed) and deformation will happen. Note bugs in the trackball mean that this behaviour isn't great, but you get the general idea. 
- Mess around with this until the application crashes!

Have fun. I hope this is vaguely useful!
 - Richard
