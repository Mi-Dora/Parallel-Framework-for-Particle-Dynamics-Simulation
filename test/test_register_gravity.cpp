#include "particle.h"
#include "user_def/gravity_particle.h"

int main(int argc, char** argv) {
    particle_t* particles;
    int nParticles;
    const std::string inputFilename = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/particle_init.txt";
    const std::string outputFilename = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/particle_iter1.txt";
    load_particles(inputFilename, &particles, &nParticles);
    
    for(int n=0; n<nParticles; n++) {
        registerUpdateAccelerationFn(particles+n, gravityUpdateAcceleration);
    }
    for(int i=0; i<nParticles; i++) {
        for(int j=0; j<nParticles; j++) {
            if(i==j)
                continue;
            updateAcceleration(particles+i, particles+j);
        }
    }
    for(int i=0; i<nParticles; i++) {
        updateVelocity(particles+i, 1e-2);
        updatePosition(particles+i, 1e-2);
    }
    store_particles(outputFilename, particles, nParticles);
    free(particles->position);
    free(particles->velocity);
    free(particles->acceleration);
    free(particles->features);
    free(particles);
    return 0;
}
