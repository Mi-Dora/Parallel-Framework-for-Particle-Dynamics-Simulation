#include "particle.h"

int main(int argc, char** argv) {
    particle_t* particles;
    int nParticles;
    std::string fLoad = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/particle_init.txt";
    std::string fStore = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/particle_1.txt";
    std::string fStore2 = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/particle_2.txt";
    std::string fStore3 = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/particle_3.txt";
    std::string fStore4 = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/particle_4.txt";
    std::string fStore5 = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/particle_5.txt";


    load_particles(fLoad, &particles, &nParticles);
    store_particles(fStore, particles, nParticles);
    free(particles->position);
    free(particles->velocity);
    free(particles->acceleration);
    free(particles->features);
    free(particles);

    load_particles(fStore, &particles, &nParticles);
    store_particles(fStore2, particles, nParticles);
    free(particles->position);
    free(particles->velocity);
    free(particles->acceleration);
    free(particles->features);
    free(particles);

    load_particles(fStore2, &particles, &nParticles);
    store_particles(fStore3, particles, nParticles);
    free(particles->position);
    free(particles->velocity);
    free(particles->acceleration);
    free(particles->features);
    free(particles);

    chunk_particles_t* chunkParticles;
    load_particles(fStore3, &chunkParticles);
    store_particles(fStore4, chunkParticles);
    free_particles(chunkParticles);

    load_particles(fStore4, &chunkParticles);
    store_particles(fStore5, chunkParticles->particles, chunkParticles->nParticle);
    free_particles(chunkParticles->particles);
    chunkParticles->particles = nullptr;
    chunkParticles->nParticle = 0;

    return 0;
}
