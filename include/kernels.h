#ifndef CUDA_KERNELS_H
#define CUDA_KERNELS_H

#include <cuda_runtime.h>
#include <cstdint>
#include "particle.h"
#include "device.h"

// void kernel(dim3& grid_size, dim3& block_size, size_t& shmem_size, 
//             double* d_pos, double* d_vel, double* d_acc, double* d_feats,
//             uint64_t* d_n_particle, uint64_t* d_n_dim, uint64_t* d_n_feat, double* d_timestep);


__global__ 
void compute_kernel(double *__restrict__ pos, 
                    double *__restrict__ vel,
                    double *__restrict__ acc,
                    double *__restrict__ feats,
                    uint64_t *__restrict__ n_particle,
                    uint64_t *__restrict__ n_dim,
                    uint64_t *__restrict__ n_feat,
                    double *__restrict__ timestep);


// __global__ 
// void scatter_kernel(particle_t const *__restrict__ particles);
#endif