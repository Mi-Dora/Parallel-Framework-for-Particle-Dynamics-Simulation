cmake_minimum_required(VERSION 3.20)
project(Parallel_Framework_for_Particle_Dynamics_Simulation)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -fopenmp -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp")

include_directories(include)

add_subdirectory(src)
add_subdirectory(test)
