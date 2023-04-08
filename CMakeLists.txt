cmake_minimum_required(VERSION 3.24)
project(Parallel_Framework_for_Particle_Dynamics_Simulation)

set(CMAKE_CXX_STANDARD 14)

include_directories(include)

add_executable(Parallel_Framework_for_Particle_Dynamics_Simulation
        include/base_particle.h
        include/kernel.h
        include/user_particle.h
        include/world.h
        src/base_particle.cpp
        src/compute_acce.cu
        src/kernel.cu
        src/main.cpp
        src/user_particle.cpp
        src/world.cpp)
