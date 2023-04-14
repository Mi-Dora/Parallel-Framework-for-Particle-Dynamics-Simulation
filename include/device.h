
#ifndef DEVICE_H
#define DEVICE_H

#include <cstdint>
#include <cuda_runtime.h>


void device_allocate_init(double** d_pos, double** d_vel, double** d_acc, double** d_feats,
                            uint64_t** d_n_particle, uint64_t** d_n_dim, uint64_t** d_n_feat, double** d_world_size, 
                            double* pos, double* vel, double* acc, double* feats, 
                            std::uint64_t n_particle, std::uint64_t n_dim, std::uint64_t n_feat, double world_size);

void kernel_params_init(dim3 &grid_size, dim3 &block_size, size_t &shmem_size, int grid_shape);

void device_free(double** d_pos, double** d_vel, double** d_acc, double** d_feats,
                            uint64_t** d_n_particle, uint64_t** d_n_dim, uint64_t** d_n_feat, double** d_world_size);

#endif