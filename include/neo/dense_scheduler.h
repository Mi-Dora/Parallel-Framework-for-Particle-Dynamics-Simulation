#ifndef __NEO_DENSE_SCHEDULER_H__
#define __NEO_DENSE_SCHEDULER_H__
#include <mpi.h>
#include <vector>
#include "neo/particle.h"

typedef struct {
    int rankX;
    int rankY;
    int gridX;
    int gridY;
    MPI_Comm reduceXComm;
    MPI_Comm reduceYComm;
    MPI_Comm scatterXComm;
    MPI_Comm scatterYComm;
} topology_t;

void setupCommunicators(topology_t* topology, int gridX, int gridY);
int minimalPaddingSize(int oriSize, int gridX, int gridY, int kernelSize);

void padding(std::vector<particle_t>* particles, int paddedSize);
void shrink(std::vector<particle_t>* particles);

void scatter(std::vector<particle_t>* localParticlesX, 
            std::vector<particle_t>* localParticlesY,
            std::vector<particle_t>* globalParticles,
            topology_t* topology);
void reduce(std::vector<particle_t>* localParticlesX, 
            std::vector<particle_t>* localParticlesY, 
            topology_t* topology);
void gather(std::vector<particle_t>* localParticlesX,
            std::vector<particle_t>* localParticlesY,
            std::vector<particle_t>* globalParticles,
            topology_t* topology);
void update(std::vector<particle_t>* localParticlesX, 
            std::vector<particle_t>* localParticlesY,
            topology_t* topology,
            double timeStep,
            void (*updateAccelerationFn)(particle_t*, particle_t*));

#endif

