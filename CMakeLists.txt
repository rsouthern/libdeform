cmake_minimum_required(VERSION 3.8)
project(libdeform_all)

# Set the installation directory. Use "cmake -DCMAKE_INSTALL_PREFIX=<install dir>" to set the installation directory where you have write permissions.

# Too lazy to type this each time
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_INSTALL_PREFIX}/bin) 

# RPATH setup which doesn't seem to work yet. You'll need to set an LD_LIBRARY_PATH to the lib to execute the test function.
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Build the child projects.
add_subdirectory(trackball)
add_subdirectory(cholmesh)
add_subdirectory(deform)
add_subdirectory(deformwidget)
