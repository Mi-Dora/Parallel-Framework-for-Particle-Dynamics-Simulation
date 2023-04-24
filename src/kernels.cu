#include "kernels.h"

__device__
inline void get_distance_vec(const double* vecA, const double* vecB, const int ndim, double* dis_vec) {
    double* a = const_cast<double*>(vecA);
    double* b = const_cast<double*>(vecB);
    for(int n=0; n<ndim; n++) 
        *(dis_vec++) = *(b++) - *(a++);
}

__device__
inline double get_length(const double* vec, const int ndim) {
    double ret = 0.0;
    double* _vec = const_cast<double*>(vec);
    for(int n=0; n<ndim; n++) {
        ret += (*_vec) * (*_vec);
        _vec++;
    }
    return sqrt(ret);
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
    uint64_t gd = gridDim.x;
    uint64_t bid = blockIdx.x;
    uint64_t tx = threadIdx.x;
    uint64_t _n_dim = *n_dim;
    uint64_t _n_feat = *n_feat;
    double _timestep = *timestep;

    const double G = 6.67e-11;
    const double pos_min = -100.;
    const double pos_max = 100.;
    
    __shared__ double share_pos_one[BLOCKSIZE][3];
    __shared__ double share_acc_one[BLOCKSIZE][3];
    __shared__ double share_pos_another[BLOCKSIZE][3];
    __shared__ double share_feat_another[BLOCKSIZE];
    uint64_t offset = bid*blockDim.x+tx;
    share_pos_one[tx][0] = *(pos+offset*_n_dim);
    share_pos_one[tx][1] = *(pos+offset*_n_dim+1);
    share_pos_one[tx][2] = *(pos+offset*_n_dim+2);
    share_acc_one[tx][0] = *(acc+offset*_n_dim);
    share_acc_one[tx][1] = *(acc+offset*_n_dim+1);
    share_acc_one[tx][2] = *(acc+offset*_n_dim+2);

    for (uint64_t i = 0; i < gd; i++){
        offset = i*blockDim.x+tx;
        share_pos_another[tx][0] = *(pos+offset*_n_dim);
        share_pos_another[tx][1] = *(pos+offset*_n_dim+1);
        share_pos_another[tx][2] = *(pos+offset*_n_dim+2);

        share_feat_another[tx] = *(feats+offset*_n_feat);
        __syncthreads();
        double dis_vec[3];
        for (uint64_t j = 0; j < BLOCKSIZE; j++){
            if (i == bid && j == tx) continue;
            get_distance_vec(share_pos_one[tx], share_pos_another[j], _n_dim, dis_vec);
            double distant = get_length(dis_vec, _n_dim);

            distant = distant * distant * distant;
            distant = G / distant * share_feat_another[j];
            for(uint64_t n=0; n < _n_dim; n++)
                share_acc_one[tx][n] += dis_vec[n] * distant;

        }
    }
    offset = bid*blockDim.x+tx;
    for (int n = 0; n < _n_dim; n++)
        *(acc+offset*_n_dim+n) = share_acc_one[tx][n];

    for (int n = 0; n < _n_dim; n++)
        *(vel+offset*_n_dim+n) += share_acc_one[tx][n]*_timestep;

    for (int n = 0; n < _n_dim; n++){
        *(pos+offset*_n_dim+n) += *(vel+offset*_n_dim+n)*_timestep;
        if (*(pos+offset*_n_dim+n) < pos_min) *(pos+offset*_n_dim+n) += pos_max-pos_min;
        if (*(pos+offset*_n_dim+n) > pos_max) *(pos+offset*_n_dim+n) -= pos_max-pos_min;
    }

}
