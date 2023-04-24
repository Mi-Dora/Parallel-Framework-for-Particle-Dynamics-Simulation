#include "particle.h"
#include "user_def/gravity_particle.h"
#include <cuda_runtime.h>
#include <chrono>
#include <cassert>

using namespace std;

int main(int argc, char* argv[]) {

    std::string n_par = argv[1];
    std::string n_iter = argv[2];
    const std::string inputFilename_cuda = "/storage/home/hcocice1/ydu340/Particle-new/data/vis/cuda_" + n_par + "_" + n_iter + ".txt";
    const std::string inputFilename_seq = "/storage/home/hcocice1/ydu340/Particle-new/data/vis/seq_" + n_par + "_" + n_iter + ".txt";

    chunk_particles_t* cuda_chunkParticles = nullptr;
    chunk_particles_t* seq_chunkParticles = nullptr;
    load_particles(inputFilename_cuda, &cuda_chunkParticles);
    load_particles(inputFilename_seq, &seq_chunkParticles);

    double *cuda_pos, *cuda_vel, *cuda_acc, *cuda_feats;
    double *seq_pos, *seq_vel, *seq_acc, *seq_feats;

    std::uint64_t n_particle, n_dim, n_feat;
    n_dim = cuda_chunkParticles->particles->ndim;
    n_particle = cuda_chunkParticles->nParticle;

    cuda_pos = cuda_chunkParticles->particles->position;
    cuda_vel = cuda_chunkParticles->particles->velocity;
    cuda_acc = cuda_chunkParticles->particles->acceleration;

    seq_pos = seq_chunkParticles->particles->position;
    seq_vel = seq_chunkParticles->particles->velocity;
    seq_acc = seq_chunkParticles->particles->acceleration;


    for (uint64_t i = 0; i < n_particle*n_dim; i++){
        assert(seq_pos[i] == cuda_pos[i]);
        assert(seq_vel[i] == cuda_vel[i]);
        assert(seq_acc[i] == cuda_acc[i]);
    }
    printf("Pass Correctness Validation!\n");
    free_particles(cuda_chunkParticles);
    free_particles(seq_chunkParticles);

    return 0;
}