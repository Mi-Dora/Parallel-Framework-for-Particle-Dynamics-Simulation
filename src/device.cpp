
#include "device.h"
#include "errs.h"
#include <cuda_runtime.h>
#include <cstdint>
#include "particle.h"

void device_allocate_init(double** d_pos, double** d_vel, double** d_acc, double** d_feats,
                            std::uint64_t** d_n_particle, std::uint64_t** d_n_dim, std::uint64_t** d_n_feat, double** d_timestep, 
                            chunk_particles_t* chunkParticles, const double timestep) {

    double *pos, *vel, *acc, *feats;

    std::uint64_t n_particle, n_dim, n_feat;
    pos = chunkParticles->particles->position;
    vel = chunkParticles->particles->velocity;
    acc = chunkParticles->particles->acceleration;
    feats = chunkParticles->particles->features;
    n_dim = chunkParticles->particles->ndim;
    n_feat = chunkParticles->particles->nfeat;
    n_particle = chunkParticles->nParticle;

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
    cuErrChk(cudaMalloc((void**)d_timestep, sizeof(double)));
    cuErrChk(cudaMemcpy(*d_timestep, &timestep, sizeof(double), cudaMemcpyHostToDevice));
}

void copy2host(double* d_pos, double* d_vel, double* d_acc, double* d_feats, chunk_particles_t* chunkParticles) {

    double *pos, *vel, *acc, *feats;
    std::uint64_t n_particle, n_dim, n_feat;
    pos = chunkParticles->particles->position;
    vel = chunkParticles->particles->velocity;
    acc = chunkParticles->particles->acceleration;
    feats = chunkParticles->particles->features;
    n_dim = chunkParticles->particles->ndim;
    n_feat = chunkParticles->particles->nfeat;
    n_particle = chunkParticles->nParticle;

    cuErrChk(cudaMemcpy(pos, d_pos, n_particle * n_dim * sizeof(double), cudaMemcpyDeviceToHost));
    cuErrChk(cudaMemcpy(vel, d_vel, n_particle * n_dim * sizeof(double), cudaMemcpyDeviceToHost));
    cuErrChk(cudaMemcpy(acc, d_acc, n_particle * n_dim * sizeof(double), cudaMemcpyDeviceToHost));
    cuErrChk(cudaMemcpy(feats, d_feats, n_particle * n_feat * sizeof(double), cudaMemcpyDeviceToHost));

}

void kernel_params_init(dim3 &grid_size, dim3 &block_size, size_t &shmem_size, int n_particle) {
    
    block_size.x = BLOCKSIZE;

    grid_size.x = n_particle / BLOCKSIZE;

    shmem_size = 0;
    // shmem_size = 16*16*3*sizeof(double);
}

void device_free(double** d_pos, double** d_vel, double** d_acc, double** d_feats,
                            uint64_t** d_n_particle, uint64_t** d_n_dim, uint64_t** d_n_feat, double** d_timestep) {
    
    cuErrChk(cudaFree(*d_pos));
    cuErrChk(cudaFree(*d_vel));
    cuErrChk(cudaFree(*d_acc));
    cuErrChk(cudaFree(*d_feats));
    cuErrChk(cudaFree(*d_n_particle));
    cuErrChk(cudaFree(*d_n_dim));
    cuErrChk(cudaFree(*d_n_feat));
    cuErrChk(cudaFree(*d_timestep));
}