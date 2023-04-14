#include "kernels.h"

__device__
inline uint64_t 3dto1d(uint64_t x_idx, uint64_t y_idx, uint64_t z_idx){
    return z_idx * gridDim.x * gridDim.y + y_idx * gridDim.x + x_idx;
}

__global__ 
void scatter_kernel(double const *__restrict__ pos, 
                    double const *__restrict__ vel,
                    double const *__restrict__ acc,
                    double const *__restrict__ feats,
                    uint64_t const *__restrict__ n_particle,
                    uint64_t const *__restrict__ n_dim,
                    uint64_t const *__restrict__ n_feat,
                    double const *__restrict__ world_size,
                    uint64_t const *__restrict__ counts
                    ){
    uint64_t bid = blockIdx.z * gridDim.x * gridDim.y + blockIdx.y * gridDim.x + blockIdx.x;
    uint64_t tid = bid * blockDim.x + threadIdx.x;
    uint64_t n_thread = gridDim.x * gridDim.y * gridDim.z * blockDim.x;
    double edge_len = world_size / gridDim.x;
    // double x_min = blockIdx.x * edge_len;
    // double y_min = blockIdx.y * edge_len;
    // double z_min = blockIdx.z * edge_len;
    // double x_max = x_min + edge_len;
    // double y_max = y_min + edge_len;
    // double z_max = z_min + edge_len;
    uint64_t x_idx, y_idx, z_idx, count;

    // count for particles in each grid
    for (int i = tid; i < n_particle; i += n_thread){
        x_idx = pos[i*n_dim] / edge_len;
        y_idx = pos[i*n_dim+1] / edge_len;
        z_idx = pos[i*n_dim+2] / edge_len;
         atomicAdd(count[3dto1d(x_idx, y_idx, z_idx)], 1);
    }

    // 
    if (blockIdx.x == 0){

    }



}

__global__ 
void compute_kernel(particle_t const *__restrict__ grided_particles,
        uint64_t const *__restrict__ offset);
