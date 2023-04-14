#include "particle.h"
#include "dense_scheduler.h"
#include "user_def/gravity_particle.h"

int main(int argc, char** argv) {
    const std::string inputFilename = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/data/data10000.txt";
    const std::string output1Filename = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/data/beforepad10000.txt";
    const std::string outputFilename = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/data/padding10000.txt";
    // const double timeStep = 0.1;

    chunk_particles_t* chunkParticles = nullptr;
    load_particles(inputFilename, &chunkParticles);
    registerUpdateAccelerationFns(chunkParticles, gravityUpdateAcceleration);

    for(int iter=0; iter<0; iter++)
        update_particles(chunkParticles, 0.1);

    printf("%lx\n", chunkParticles->particles);

    store_particles(output1Filename, chunkParticles);
    padding(chunkParticles, 3);

    for(int iter=0; iter<0; iter++)
        update_particles(chunkParticles, 0.1);

    printf("%lx\n", chunkParticles->particles);
    particle_t* test = chunkParticles->particles + 100;
    printf("%lx\n", test);
    printf("%lx %lx %lx %lx\n", test->position, test->velocity, test->acceleration, test->features);
    std::printf("nParticles=%d\n", chunkParticles->nParticle);
    std::printf("nparticles=%d id=%ld enabled=%d pos=%f vel=%f acc=%f feat=%f\n", chunkParticles->nParticle, test->id, test->enabled, *(test->position), *(test->velocity), *(test->acceleration), *(test->features));
    shrink(chunkParticles);
    
    store_particles(outputFilename, chunkParticles);
    free_particles(chunkParticles);
    free(chunkParticles);
    chunkParticles = nullptr;
    return 0;
}