#include "kernels.h"
#include "device.h"
#include "particle.h"
#include <cuda_runtime.h>



int main(int argc, char** argv) {

    const std::string inputFilename = "/storage/home/hcocice1/ydu340/Particle-new/data/data16384.txt";
    const double timeStep = 0.1;

    chunk_particles_t* chunkParticles = nullptr;
    load_particles(inputFilename, &chunkParticles);


    double *d_pos, *d_vel, *d_acc, *d_feats, *d_timestep;
    uint64_t *d_n_particle, *d_n_dim, *d_n_feat;
    device_allocate_init(&d_pos, &d_vel, &d_acc, &d_feats, &d_n_particle, &d_n_dim, &d_n_feat, &d_timestep, chunkParticles, timeStep);

    dim3 grid_size, block_size;
    size_t shmem_size;
    kernel_params_init(grid_size, block_size, shmem_size, chunkParticles->nParticle);
    
    cudaEvent_t start, stop;
    cuErrChk(cudaEventCreate(&start));
    cuErrChk(cudaEventCreate(&stop));

    float cur_time = 0.;
    float tt_time = 0;
    for(int iter=0; iter<1001; iter++) {
        
        cuErrChk(cudaEventRecord(start, NULL));
        compute_kernel <<< grid_size, block_size, shmem_size >>> (
                d_pos, d_vel, d_acc, d_feats, d_n_particle, d_n_dim, d_n_feat, d_timestep);
        cuErrChk(cudaEventRecord(stop, NULL));
        cuErrChk(cudaEventSynchronize(stop));
        cur_time = 0;
        cuErrChk(cudaEventElapsedTime(&cur_time, start, stop));
        tt_time += cur_time;
        if(iter%50 == 0) {
            printf("iter=%d, n_particle=%d, avg compute time per iter=%.6f ms\n", iter, chunkParticles->nParticle, tt_time/(iter+1));
            copy2host(d_pos, d_vel, d_acc, d_feats, chunkParticles);
            std::string outputFilename = "/storage/home/hcocice1/ydu340/Particle-new/data/16384_"  + std::to_string(iter) + ".txt";
            store_particles(outputFilename, chunkParticles);
        }
    }

    free_particles(chunkParticles);
    return 0;
}