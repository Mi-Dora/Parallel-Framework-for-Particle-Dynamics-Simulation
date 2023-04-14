#ifndef CUDA_KERNELS_H
#define CUDA_KERNELS_H

#include <cuda_runtime.h>
#include <cstdint>
#include "particle.h"

__global__ 
void scatter_kernel(particle_t const *__restrict__ particles);

__global__ 
void compute_kernel(particle_t const *__restrict__ grided_particles,
        uint64_t const *__restrict__ offset);
#endif