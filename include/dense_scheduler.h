#ifndef __DENSE_SCHEDULER_H__
#define __DENSE_SCHEDULER_H__
#include <mpi.h>
#include "particle.h"

typedef struct {
    int rankX;
    int rankY;
    int gridX;
    int gridY;
    MPI_Comm gatherXComm;
    MPI_Comm gatherYComm;
    MPI_Comm scatterXComm;
    MPI_Comm scatterYComm;
} topology_t;

void setupCommunicators(topology_t* topology, int gridX, int gridY);
void padding(particle_t** particles, int* nParticles, int minimalGrid);
void padding(chunk_particles_t* chunkParticles, int minimalGrid);

void shrink(particle_t** particles, int* nParticles);
void shrink(chunk_particles_t* chunkParticles);

#endif
