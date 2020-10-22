cmake_minimum_required(VERSION 3.8)

# Name of the project with Language
project(cholmesh)

set(CMAKE_CXX_STANDARD 11)

include_directories(include ${CMAKE_SOURCE_DIR}/include)

#find_library(CHOLMOD cholmod)
#find_library(OPENMESH OpenMesh)

# build the executable
add_library(
    cholmesh
    SHARED
    src/mesh.cpp
    src/dual.cpp
    src/cholmod_shared.cpp
)

target_link_libraries(cholmesh -lOpenMeshCore -lOpenMeshTools -lcholmod)


