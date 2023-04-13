#include "particle.h"

int main(int argc, char** argv) {
    particle_t* particles;
    int nParticles;
    std::string fLoad = "particle_init.txt";
    std::string fStore = "particle_store.txt";

    load_particles(fLoad, &particles, &nParticles);
    store_particles(fStore, particles, nParticles);
    free(particles->position);
    free(particles->velocity);
    free(particles->acceleration);
    free(particles->features);
    free(particles);

    return 0;
}
