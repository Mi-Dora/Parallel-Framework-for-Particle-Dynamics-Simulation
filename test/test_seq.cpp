#include "particle.h"
#include "user_def/gravity_particle.h"
#include <cuda_runtime.h>
#include <chrono>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {

    std::string n_par = argv[1];
    int n_iter = std::atoi(argv[2]);
    int save_interval = std::atoi(argv[3]);
    const std::string inputFilename = "/storage/home/hcocice1/ydu340/Particle-new/data/data" + n_par + ".txt";
    const double timeStep = 0.1;

    chunk_particles_t* chunkParticles = nullptr;
    load_particles(inputFilename, &chunkParticles);
    double pos_min = -100.;
    double pos_max = 100.;

    double *pos, *vel, *acc, *feats;

    std::uint64_t n_particle, n_dim, n_feat;
    pos = chunkParticles->particles->position;
    vel = chunkParticles->particles->velocity;
    acc = chunkParticles->particles->acceleration;
    feats = chunkParticles->particles->features;
    n_dim = chunkParticles->particles->ndim;
    n_feat = chunkParticles->particles->nfeat;
    n_particle = chunkParticles->nParticle;

    float cur_time = 0.;
    
    for(int iter=0; iter<n_iter; iter++) {
        printf("Seq iter %d\n", iter);
        particle_t *one = static_cast<particle_t*>(malloc(sizeof(particle_t)));
        particle_t *another = static_cast<particle_t*>(malloc(sizeof(particle_t)));

        auto start = std::chrono::steady_clock::now();

    
        one->ndim = n_dim;
        another->ndim = n_dim;
        one->nfeat = n_feat;
        another->nfeat = n_feat;
        // double 

        for (uint64_t i = 0; i < n_particle; i++){
            one->position = pos+i*n_dim;
            one->velocity = vel+i*n_dim;
            one->acceleration = acc+i*n_dim;
            one->features = feats+i*n_feat;
            for (uint64_t j = 0; j < n_particle; j++){
                if (i == j) continue;
                another->position = pos+j*n_dim;
                another->velocity = vel+j*n_dim;
                another->acceleration = acc+j*n_dim;
                another->features = feats+j*n_feat;
                gravityUpdateAcceleration(one, another);
            }
        }
        for (uint64_t i = 0; i < n_particle; i++){
            one->position = pos+i*n_dim;
            one->velocity = vel+i*n_dim;
            one->acceleration = acc+i*n_dim;
            for(uint64_t n = 0; n < n_dim; n++){
                *(one->velocity+n) += *(one->acceleration+n) * timeStep;
                *(one->position+n) += *(one->velocity+n) * timeStep;
                if (*(one->position+n) < pos_min) *(one->position+n) += pos_max-pos_min;
                if (*(one->position+n) > pos_max) *(one->position+n) -= pos_max-pos_min;
            }
        }
        auto end = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        if(iter%save_interval == 0) {
            cout << "Sequential: n_particle = " << n_particle << ", compute time per iter = " << elapsed.count() << "ms" << endl;
            // printf("iter=%d, n_particle=%d, avg compute time per iter=%.6f ms\n", iter, chunkParticles->nParticle, tt_time/(iter+1));
            std::string outputFilename = "/storage/home/hcocice1/ydu340/Particle-new/data/vis/seq_" + n_par + "_"  + std::to_string(iter) + ".txt";
            store_particles(outputFilename, chunkParticles);
        }
    }

    free_particles(chunkParticles);
    return 0;
}