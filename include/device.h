
#ifndef DEVICE_H
#define DEVICE_H

#include <cstdint>
#include <cuda_runtime.h>
#include "errs.h"
#include "particle.h"

#define BLOCKSIZE 128

void device_allocate_init(double** d_pos, double** d_vel, double** d_acc, double** d_feats,
                            std::uint64_t** d_n_particle, std::uint64_t** d_n_dim, std::uint64_t** d_n_feat, double** d_timestep, 
                            chunk_particles_t* chunkParticles, const double timestep);

void copy2host(double* d_pos, double* d_vel, double* d_acc, double* d_feats, chunk_particles_t* chunkParticles);

void kernel_params_init(dim3 &grid_size, dim3 &block_size, size_t &shmem_size, int n_particle);

void device_free(double** d_pos, double** d_vel, double** d_acc, double** d_feats,
                            uint64_t** d_n_particle, uint64_t** d_n_dim, uint64_t** d_n_feat, double** d_timestep);

#endif