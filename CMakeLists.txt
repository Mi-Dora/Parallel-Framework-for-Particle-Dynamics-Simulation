cmake_minimum_required(VERSION 3.20)
project(Parallel_Framework_for_Particle_Dynamics_Simulation)

set(CMAKE_CXX_STANDARD 14)

include_directories(include)

add_executable(Parallel_Framework_for_Particle_Dynamics_Simulation
        src/particle.cpp 
        src/utils.cpp 
        src/world.cpp 
        src/user_def/gravity_particle.cpp 
        src/kernels.cu
)
