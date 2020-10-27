cmake_minimum_required(VERSION 3.8)

# Name of the project with Language
project(cholmesh)

set(CMAKE_CXX_STANDARD 11)

include_directories(include ${CMAKE_SOURCE_DIR}/include $ENV{OPENMESH_PATH}/include $ENV{SUITESPARSE_PATH}/include $ENV{SUITESPARSE_PATH})

find_library(OPENMESH_LIBRARIES OpenMeshCore OpenMeshTools PATHS $ENV{OPENMESH_PATH}/lib)
find_library(SUITESPARSE_LIBRARIES cholmod PATHS $ENV{SUITESPARSE_PATH}/lib)

# build the executable
add_library(
    cholmesh
    SHARED
    src/mesh.cpp
    src/dual.cpp
    src/cholmod_shared.cpp
)

target_link_libraries(cholmesh ${OPENMESH_LIBRARIES} ${SUITESPARSE_LIBRARIES})


