#ifndef __NEO_CUTOFF_SCHEDULER_H__
#define __NEO_CUTOFF_SCHEDULER_H__
#include "neo/particle.h"
#include <vector>

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


void divideGroupN(std::vector<particle_t>* allParticles, 
                std::vector<std::vector<particle_t>>* groupParticles,
                int ndim,
                int ngroup);

void dispatch(std::vector<particle_t>* allParticles, 
                std::vector<std::vector<particle_t>>* localParticleGroups,
                topology_t* topology
                );

#endif
