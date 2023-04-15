#ifndef CUDA_KERNELS_H
#define CUDA_KERNELS_H

#include <cuda_runtime.h>
#include <cstdint>
#include "particle.h"


__global__ 
void compute_kernel(double const *__restrict__ pos, 
                    double const *__restrict__ vel,
                    double const *__restrict__ acc,
                    double const *__restrict__ feats,
                    uint64_t const *__restrict__ n_particle,
                    uint64_t const *__restrict__ n_dim,
                    uint64_t const *__restrict__ n_feat,
                    double const *__restrict__ timestep);
                    
// __global__ 
// void scatter_kernel(particle_t const *__restrict__ particles);
#endif