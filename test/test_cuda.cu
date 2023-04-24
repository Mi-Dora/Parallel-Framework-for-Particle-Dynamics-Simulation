#include "kernels.h"
#include "device.h"
#include "particle.h"
#include <cuda_runtime.h>
using namespace std;



int main(int argc, char** argv) {

    std::string n_par = argv[1];
    int n_iter = std::atoi(argv[2]);
    int save_interval = std::atoi(argv[3]);
    const std::string inputFilename = "/storage/home/hcocice1/ydu340/Particle-new/data/data" + n_par + ".txt";
    const double timeStep = 0.1;

    chunk_particles_t* chunkParticles = nullptr;
    load_particles(inputFilename, &chunkParticles);


    double *d_pos, *d_vel, *d_acc, *d_feats, *d_timestep;
    uint64_t *d_n_particle, *d_n_dim, *d_n_feat;
    float cur_time = 0.;
    cudaEvent_t start, stop;
    cuErrChk(cudaEventCreate(&start));
    cuErrChk(cudaEventCreate(&stop));

    cuErrChk(cudaEventRecord(start, NULL));
    clock_t start_c = clock();
    device_allocate_init(&d_pos, &d_vel, &d_acc, &d_feats, &d_n_particle, &d_n_dim, &d_n_feat, &d_timestep, chunkParticles, timeStep);
    clock_t end_c   = clock();
    
    cuErrChk(cudaEventRecord(stop, NULL));
    cuErrChk(cudaEventSynchronize(stop));
    cur_time = 0;
    cuErrChk(cudaEventElapsedTime(&cur_time, start, stop));

    printf("n_particle=%d, allocate memory time = %.6f ms\n", chunkParticles->nParticle, cur_time);
    dim3 grid_size, block_size;
    size_t shmem_size;
    kernel_params_init(grid_size, block_size, shmem_size, chunkParticles->nParticle);
    
    

    
    float tt_time = 0.;
    // float mem_time = 0.;
    for(int iter=0; iter<n_iter; iter++) {
        
        cuErrChk(cudaEventRecord(start, NULL));
        compute_kernel <<< grid_size, block_size, shmem_size >>> (
                d_pos, d_vel, d_acc, d_feats, d_n_particle, d_n_dim, d_n_feat, d_timestep);
        cuErrChk(cudaEventRecord(stop, NULL));
        cuErrChk(cudaEventSynchronize(stop));
        cur_time = 0.;
        cuErrChk(cudaEventElapsedTime(&cur_time, start, stop));
        tt_time += cur_time;
        if(iter%save_interval == 0) {
            printf("iter=%d, n_particle=%d, avg compute time per iter=%.6f ms\n", iter, chunkParticles->nParticle, tt_time/(iter+1));
            cuErrChk(cudaEventRecord(start, NULL));

            copy2host(d_pos, d_vel, d_acc, d_feats, chunkParticles);

            cuErrChk(cudaEventRecord(stop, NULL));
            cuErrChk(cudaEventSynchronize(stop));
            cur_time = 0;
            cuErrChk(cudaEventElapsedTime(&cur_time, start, stop));
            printf("n_particle=%d, copy back memory time = %.6f ms\n", chunkParticles->nParticle, cur_time);
            std::string outputFilename = "/storage/home/hcocice1/ydu340/Particle-new/data/cuda_" + n_par + "_"   + std::to_string(iter) + ".txt";
            store_particles(outputFilename, chunkParticles);
        }
    }

    free_particles(chunkParticles);
    return 0;
}