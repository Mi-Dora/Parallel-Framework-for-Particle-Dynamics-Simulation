#include "kernels.h"


int main(int argc, char** argv) {
    return 0;
    // const std::string inputFilename = "../data/data10000.txt";
    // const double timeStep = 0.1;
    // const double cell_size = 10;

    // chunk_particles_t* chunkParticles = nullptr;
    // load_particles(inputFilename, &chunkParticles);

    // double *pos, *vel, *acc, *feats;
    // int n_particle, n_dim, n_feat;
    // pos = chunkParticles->particles->position;
    // vel = chunkParticles->particles->velocity;
    // acc = chunkParticles->particles->acceleration;
    // feats = chunkParticles->particles->features;
    // n_dim = chunkParticles->particles->ndim;
    // n_feat = chunkParticles->particles->nfeat;
    // n_particle = chunkParticles->nParticles;

    // double *d_pos, *d_vel, *d_acc, *d_feats, *d_timestep;
    // uint64_t *d_n_particle, *d_n_dim, *d_n_feat, 
    // device_allocate_init(&d_pos, &d_vel, &d_acc, &d_feats, 
    // &d_n_particle, &d_n_dim, &d_n_feat, &d_timestep, 
    // pos, vel, acc, feats, n_particle, n_dim, n_feat, timeStep);

    // dim3 grid_size, block_size;
    // size_t shmem_size;
    // kernel_params_init(grid_size, block_size, shmem_size, cell_size);
    
    // cudaEvent_t start, stop;
    // cuErrChk(cudaEventCreate(&start));
    // cuErrChk(cudaEventCreate(&stop));
    // cublasHandle_t handle;
    // cuErrChk(cublasCreate(&handle));

    // float cur_time = 0.;
    // for(int iter=0; iter<1000; iter++) {
    //     printf("iter=%d\n", iter);
    //     cuErrChk(cudaEventRecord(start, NULL));
    //     compute_kernel <<< grid_size, block_size, shmem_size >>> (
    //             d_pos, d_vel, d_acc, d_feats, d_n_particle, d_n_dim, d_n_feat, d_timestep);
    //     cuErrChk(cudaEventRecord(stop, NULL));
    //     cuErrChk(cudaEventSynchronize(stop));
    //     cur_time = 0;
    //     cuErrChk(cudaEventElapsedTime(&cur_time, start, stop));

    //     if(iter%50 == 0) {
    //         std::string outputFilename = "../data/test_updates/"  + std::to_string(iter) + ".txt";
    //         store_particles(outputFilename, chunkParticles);
    //     }
    // }

    // free_particles(chunkParticles);
    
}