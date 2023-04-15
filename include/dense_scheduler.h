#ifndef __DENSE_SCHEDULER_H__
#define __DENSE_SCHEDULER_H__
#include <mpi.h>
#include "particle.h"

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
void padding(particle_t** particles, int* nParticles, int minimalGrid);
void padding(chunk_particles_t* chunkParticles, int minimalGrid);

void shrink(particle_t** particles, int* nParticles);
void shrink(chunk_particles_t* chunkParticles);

void scatter(particle_t** localParticlesX, int* localNParticlesX,
            particle_t** localParticlesY, int* localNParticlesY,
            particle_t* wholeParticles, int wholeNParticles, 
            topology_t* topology);
void scatter(chunk_particles_t** localChunkParticlesX, 
                chunk_particles_t** localChunkParticlesY,
                chunk_particles_t* wholeChunkParticles,
                topology_t* topology);

void reduce(particle_t* localParticlesX, int localNParticlesX, 
        particle_t* localParticlesY, int localNParticlesY, topology_t* topology);
void reduce(chunk_particles_t* localChunkParticlesX, chunk_particles_t* localChunkParticlesY, topology_t* topology);


void gather(particle_t* localParticlesX, int localNParticlesX, 
                particle_t* localParticlesY, int localNParticlesY,
                particle_t** globalParticles, int* globalNParticles, topology_t* topology);
void gather(chunk_particles_t* localChunkParticlesX, 
            chunk_particles_t* localChunkParticlesY, 
            chunk_particles_t* globalChunkParticles, 
            topology_t* topology);

#endif
