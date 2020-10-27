cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(cmake_and_cuda LANGUAGES CXX CUDA)

# Note that the CUDA_SAMPLES_PATH needs to be set as an environment variable
include_directories(include ${CMAKE_SOURCE_DIR}/include $ENV{SVDLIB_PATH} ../cholmesh/include $ENV{CUDA_SAMPLES_PATH}/Common $ENV{CUDA_DEPRECATED} $ENV{SUITESPARSE_PATH}/include $ENV{SVDLIB_PATH} $ENV{OPENMESH_PATH}/include)

# Uncomment the line below to get verbose compilation information
#set(CMAKE_VERBOSE_MAKEFILE ON)

find_package(CUDA)
find_library(OPENMESH_LIBRARIES OpenMeshCore OpenMeshTools PATHS $ENV{OPENMESH_PATH}/lib)
find_library(SUITESPARSE_LIBRARIES cholmod PATHS $ENV{SUITESPARSE_PATH}/lib)

add_library(deform SHARED
src/cgaltools.cpp
src/laplacian_DUAL2.cpp
src/mattools.cpp
src/cholmod_shared.cpp
src/dual.cpp           
src/laplacian_DUAL2_cuda.cu  
src/mesh.cpp
src/cuda_helper.cu      
src/laplacian_COT.cpp  
src/laplacian_DUAL.cpp       
src/reduce.cu
src/cusp_tools.cu       
src/laplacian.cpp      
src/laplacian_DUAL_cuda.cu
)

# Request that deform is built with -std=c++11
target_compile_features(deform PUBLIC cxx_std_11)

# We need to explicitly state that we need all CUDA files in the 
# library to be built with -dc as the member functions could be called by
# other libraries and executables
set_target_properties(deform PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# Link in CUDA and CURAND
target_link_libraries(deform ${CUDA_LIBRARIES} ${CUDA_curand_LIBRARY} ${SUITESPARSE_LIBRARIES} ${OPENMESH_LIBRARIES} $ENV{SVDLIB_PATH}/libsvd.a)

