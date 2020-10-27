cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project(deformwidget LANGUAGES CXX CUDA)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

set (QTVIEWER_PATH ${CMAKE_CURRENT_SOURCE_DIR}/package/OpenMesh/Apps/QtViewer)

# include (ACGCommon)

# Comment out if you don't want all the gory detail
set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED OpenGL REQUIRED)

#find_library(CHOLMESH_LIBRARY cholmesh PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#find_library(DEFORM_LIBRARY deform PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
#find_library(TRACKBALL_LIBRARY trackball PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

# Set the OpenGL preference. Alternative option is "LEGACY"
set(OpenGL_GL_PREFERENCE GLVND)
find_package(OpenGL REQUIRED)

# find_package(GLEW REQUIRED)
find_path(GLEW_INCLUDE_DIR GL/glew.h)
find_library(GLEW_LIBRARY NAMES GLEW glew32 glew glew32s PATH_SUFFIXES lib64)

find_package(GLUT REQUIRED)

# The stuff below is so that CMake can identify QObjects and use MOC on them
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

include_directories (
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${CMAKE_CURRENT_SOURCE_DIR}/package
	${CMAKE_SOURCE_DIR}/cholmesh/include 
	${CMAKE_SOURCE_DIR}/deform/include 
	${CMAKE_SOURCE_DIR}/trackball/include 
	${CMAKE_SOURCE_DIR}/include
)

add_executable(deformwidget
	${QTVIEWER_PATH}/QGLViewerWidget.cc
	${QTVIEWER_PATH}/MeshViewerWidget.cc
	${QTVIEWER_PATH}/MeshViewerWidgetT_impl.hh
	${QTVIEWER_PATH}/MeshViewerWidget.hh
	${QTVIEWER_PATH}/QGLViewerWidget.hh
	${QTVIEWER_PATH}/MeshViewerWidgetT.hh
	src/choltrackball.cpp 
	src/deformwidget.cpp 
	src/framebuffer.cpp 
	src/glController.cpp 
	src/main.cpp 
	src/vertexcolormap.cpp
	include/choltrackball.h
	include/deformwidget.h
	include/framebuffer.h
	include/glController.h
	include/vertexcolormap.h
)

target_link_libraries (deformwidget
  OpenMeshCore
  OpenMeshTools
  cholmod
  ${GLEW_LIBRARIES}
  ${GLUT_LIBRARIES}
  ${OPENGL_LIBRARIES}
  #${CHOLMESH_LIBRARY}
  #${DEFORM_LIBRARY}
  #${TRACKBALL_LIBRARY}
  cholmesh
  deform
  trackball
  Qt5::Widgets
  Qt5::OpenGL
  gmp
)

