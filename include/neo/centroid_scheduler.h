#ifndef __NEO_CENTROID_SCHEDULER_H__
#define __NEO_CENTROID_SCHEDULER_H__
#include "neo/particle.h"
#include <vector>
#include <mpi.h>

typedef struct {
    int rankX;
    int rankY;
    int rankZ;
    int gridX;
    int gridY;
    int gridZ;
    int subgridX;
    int subgridY;
    int subgridZ;
    MPI_Comm xEdgeComm;
    MPI_Comm yEdgeComm;
    MPI_Comm zEdgeComm;
} topology_t;

void setupCommunicators(topology_t* topology, int gridX, int gridY, int gridZ, int subGridX, int subGridY, int subGridZ);



void divideGroupN(std::vector<particle_t>* allParticles, 
                std::vector<std::vector<particle_t>>* groupParticles,
                int ndim,
                int ngroup);

void dispatch(std::vector<particle_t>* allParticles, 
                std::vector<std::vector<particle_t>>* localParticleGroups,
                topology_t* topology
                );

void gather(std::vector<std::vector<particle_t>>* localParticleGroups, 
            std::vector<particle_t>* allParticles, 
            topology_t* topology);

void getCentroid(std::vector<particle_t>* particles, particle_t** centroid, int weightFeat);

#endif
