cmake_minimum_required(VERSION 3.20)
project(Parallel_Framework_for_Particle_Dynamics_Simulation)

include(CheckLanguage)
check_language(CUDA)
enable_language(CUDA)

option(ENABLE_GPU "Enable GPU acceleration for parallel particle simulation" ON)

if(ENABLE_GPU)
    add_definitions(-DENABLE_GPU)
endif()

# set(CMAKE_C_COMPILER gcc)
# set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -fopenmp -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -fopenmp")

set(CMAKE_CUDA_STANDARD 11)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
set(CMAKE_CUDA_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CUDA_ARCHITECTURES 35 50 70 72)
set(CMAKE_CUDA_FLAGS_DEBUG "-G")
set(CMAKE_CUDA_FLAG_RELEASE "-Wno-deprecated-gpu-targets")
set(CMAKE_CUDA_ARCHITECTURES 70 72)

find_package(CUDAToolkit REQUIRED)

include_directories(include)

add_subdirectory(src)
add_subdirectory(test)
