
#include <device.h>
#include <errs.h>
#include <cuda_runtime.h>
#include <cstdint>
#include "particle.h"

void device_allocate_init(double** d_pos, double** d_vel, double** d_acc, double** d_feats,
                            uint64_t** d_n_particle, uint64_t** d_n_dim, uint64_t** d_n_feat, double** d_world_size, 
                            double* pos, double* vel, double* acc, double* feats, 
                            std::uint64_t n_particle, std::uint64_t n_dim, std::uint64_t n_feat, double world_size) {

    cuErrChk(cudaMalloc((void**)d_pos, n_particle* n_dim * sizeof(double)));
    cuErrChk(cudaMemcpy(*d_pos, pos, n_particle * n_dim * sizeof(double), cudaMemcpyHostToDevice));
    cuErrChk(cudaMalloc((void**)d_vel, n_particle* n_dim * sizeof(double)));
    cuErrChk(cudaMemcpy(*d_vel, vel, n_particle * n_dim * sizeof(double), cudaMemcpyHostToDevice));
    cuErrChk(cudaMalloc((void**)d_acc, n_particle* n_dim * sizeof(double)));
    cuErrChk(cudaMemcpy(*d_acc, acc, n_particle * n_dim * sizeof(double), cudaMemcpyHostToDevice));
    cuErrChk(cudaMalloc((void**)d_feats, n_particle* n_feat * sizeof(double)));
    cuErrChk(cudaMemcpy(*d_feats, feats, n_particle * n_feat * sizeof(double), cudaMemcpyHostToDevice));

    cuErrChk(cudaMalloc((void**)d_n_particle, sizeof(std::uint64_t)));
    cuErrChk(cudaMemcpy(*d_n_particle, &n_particle, sizeof(std::uint64_t), cudaMemcpyHostToDevice));
    cuErrChk(cudaMalloc((void**)d_n_dim, sizeof(std::uint64_t)));
    cuErrChk(cudaMemcpy(*d_n_dim, &n_dim, sizeof(std::uint64_t), cudaMemcpyHostToDevice));
    cuErrChk(cudaMalloc((void**)d_n_feat, sizeof(std::uint64_t)));
    cuErrChk(cudaMemcpy(*d_n_feat, &n_feat, sizeof(std::uint64_t), cudaMemcpyHostToDevice));
    cuErrChk(cudaMalloc((void**)d_world_size, sizeof(double)));
    cuErrChk(cudaMemcpy(*d_world_size, &world_size, sizeof(double), cudaMemcpyHostToDevice));
}

void kernel_params_init(dim3 &grid_size, dim3 &block_size, size_t &shmem_size, int cell_size) {
    
    block_size.x = 1024;
    // block_size.y = 32;

    grid_size.x = grid_shape;
    grid_size.y = grid_shape;
    grid_size.z = grid_shape;

    shmem_size = 0;
}

void device_free(double** d_pos, double** d_vel, double** d_acc, double** d_feats,
                            uint64_t** d_n_particle, uint64_t** d_n_dim, uint64_t** d_n_feat, double** d_world_size) {
    
    cuErrChk(cudaFree(*d_pos));
    cuErrChk(cudaFree(*d_vel));
    cuErrChk(cudaFree(*d_acc));
    cuErrChk(cudaFree(*d_feats));
    cuErrChk(cudaFree(*d_n_particle));
    cuErrChk(cudaFree(*d_n_dim));
    cuErrChk(cudaFree(*d_n_feat));
    cuErrChk(cudaFree(*d_world_size));
}