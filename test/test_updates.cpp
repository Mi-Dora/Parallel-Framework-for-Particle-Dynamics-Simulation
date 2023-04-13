#include "particle.h"
#include "user_def/gravity_particle.h"

int main(int argc, char** argv) {
    const std::string inputFilename = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/particle_init.txt";
    const double timeStep = 0.1;

    chunk_particles_t* chunkParticles = nullptr;
    load_particles(inputFilename, &chunkParticles);
    registerUpdateAccelerationFns(chunkParticles, gravityUpdateAcceleration);

    for(int iter=0; iter<1000; iter++) {
        update_particles(chunkParticles, timeStep);
        if(iter%50 == 0) {
            std::string outputFilename = inputFilename + "." + std::to_string(iter);
            store_particles(outputFilename, chunkParticles);
        }
    }
    free_particles(chunkParticles);
    
}