cmake_minimum_required(VERSION 3.8)

# Name of the project with Language
project(trackball)

set(CMAKE_CXX_STANDARD 11)

include_directories(include ${CMAKE_SOURCE_DIR}/include)

# build the executable
add_library(trackball
    SHARED
    src/glTrackball.cpp 
    src/gmMat3.cpp 
    src/gmMat4.cpp 
    src/gmQuat.cpp 
    src/Trackball.cpp
)

# target_link_libraries(cholmesh -lOpenMeshCore -lOpenMeshTools -lcholmod)


