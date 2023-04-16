#include "kernels.h"


// __device__
// inline uint64_t 3dto1d(uint64_t x_idx, uint64_t y_idx, uint64_t z_idx){
//     return z_idx * gridDim.x * gridDim.y + y_idx * gridDim.x + x_idx;
// }

__device__
void get_distance_vec(const double* vecA, const double* vecB, const int ndim, double** vecRet) {
    double* ret;
    double* a = const_cast<double*>(vecA);
    double* b = const_cast<double*>(vecB);
    ret = static_cast<double*>(malloc(sizeof(double)*ndim));
    *vecRet = ret;
    for(int n=0; n<ndim; n++) 
        *(ret++) = *(a++) - *(b++);
}

__device__
double get_length(const double* vec, const int ndim) {
    double ret = 0.0;
    double* _vec = const_cast<double*>(vec);
    for(int n=0; n<ndim; n++) {
        ret += (*_vec) * (*_vec);
        _vec++;
    }
    return sqrt(ret);
}

__device__
void d_gravityUpdateAcceleration(particle_t* one, particle_t* another) {
    const double G = 6.67e-11;
    // double* mass1 = one->features;
    double* mass2 = another->features;

    double* positionDiff = nullptr;
    // printf("%d ", one->ndim);
    get_distance_vec(one->position, another->position, one->ndim, &positionDiff);
    double distant = get_length(positionDiff, one->ndim);
    distant = distant * distant * distant;
    distant = G / distant * (*mass2);
    double* acc = one->acceleration;
    double* positionDiff_ = positionDiff;
    for(int n=0; n<one->ndim; n++)
        *(acc++) += *(positionDiff_++) * distant;
    free(positionDiff);
}

__global__ 
void compute_kernel(double *__restrict__ pos, 
                    double *__restrict__ vel,
                    double *__restrict__ acc,
                    double *__restrict__ feats,
                    uint64_t *__restrict__ n_particle,
                    uint64_t *__restrict__ n_dim,
                    uint64_t *__restrict__ n_feat,
                    double *__restrict__ timestep)
{
    uint64_t bid = blockIdx.z * gridDim.x * gridDim.y + blockIdx.y * gridDim.x + blockIdx.x;
    uint64_t tid = bid * blockDim.x + threadIdx.x;
    uint64_t n_thread = gridDim.x * gridDim.y * gridDim.z * blockDim.x;
    uint64_t _n_dim = *n_dim;
    particle_t *one = static_cast<particle_t*>(malloc(sizeof(particle_t)));
    particle_t *another = static_cast<particle_t*>(malloc(sizeof(particle_t)));
    one->ndim = _n_dim;
    another->ndim = _n_dim;

    for (uint64_t i = tid; i < (*n_particle); i += n_thread){
        one->position = pos+i*_n_dim;
        one->velocity = vel+i*_n_dim;
        one->acceleration = acc+i*_n_dim;
        one->features = feats+i*(*n_feat);
        for (uint64_t j = 0; j < (*n_particle); j++){
            if (i == j) continue;
            another->position = pos+j*_n_dim;
            another->velocity = vel+j*_n_dim;
            another->acceleration = acc+j*_n_dim;
            another->features = feats+j*(*n_feat);
            d_gravityUpdateAcceleration(one, another);
        }
    }
    for (uint64_t i = tid; i < (*n_particle); i += n_thread){
        one->position = pos+i*_n_dim;
        one->velocity = vel+i*_n_dim;
        one->acceleration = acc+i*_n_dim;
        for(int n = 0; n < _n_dim; n++){
            *(one->velocity+n) += *(one->acceleration+n) * (*timestep);
            *(one->position+n) += *(one->velocity+n) * (*timestep);
        }
    }

}

// void kernel(dim3& grid_size, dim3& block_size, size_t& shmem_size, 
//             double* d_pos, double* d_vel, double* d_acc, double* d_feats,
//             uint64_t* d_n_particle, uint64_t* d_n_dim, uint64_t* d_n_feat, double* d_timestep){
//     compute_kernel <<< grid_size, block_size, shmem_size >>> (
//                 d_pos, d_vel, d_acc, d_feats, d_n_particle, d_n_dim, d_n_feat, d_timestep);
// }

// __global__ 
// void scatter_kernel(double const *__restrict__ pos, 
//                     double const *__restrict__ vel,
//                     double const *__restrict__ acc,
//                     double const *__restrict__ feats,
//                     uint64_t const *__restrict__ n_particle,
//                     uint64_t const *__restrict__ n_dim,
//                     uint64_t const *__restrict__ n_feat,
//                     double const *__restrict__ world_size,
//                     uint64_t const *__restrict__ counts
//                     ){
//     uint64_t bid = blockIdx.z * gridDim.x * gridDim.y + blockIdx.y * gridDim.x + blockIdx.x;
//     uint64_t tid = bid * blockDim.x + threadIdx.x;
//     uint64_t n_thread = gridDim.x * gridDim.y * gridDim.z * blockDim.x;
//     double edge_len = world_size / gridDim.x;
//     // double x_min = blockIdx.x * edge_len;
//     // double y_min = blockIdx.y * edge_len;
//     // double z_min = blockIdx.z * edge_len;
//     // double x_max = x_min + edge_len;
//     // double y_max = y_min + edge_len;
//     // double z_max = z_min + edge_len;
//     uint64_t x_idx, y_idx, z_idx, count;

//     // count for particles in each grid
//     for (int i = tid; i < n_particle; i += n_thread){
//         x_idx = pos[i*n_dim] / edge_len;
//         y_idx = pos[i*n_dim+1] / edge_len;
//         z_idx = pos[i*n_dim+2] / edge_len;
//          atomicAdd(count[3dto1d(x_idx, y_idx, z_idx)], 1);
//     }

//     // 
//     if (blockIdx.x == 0){

//     }



// }

// __global__ 
// void compute_cutoff_kernel(particle_t const *__restrict__ grided_particles,
//         uint64_t const *__restrict__ offset);
