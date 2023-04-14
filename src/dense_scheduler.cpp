#include "dense_scheduler.h"
#include <mpi.h>
#include <stdlib.h>
#include <cassert>

void setupCommunicators(topology_t* topology, int gridX, int gridY) {
    topology->gridX = gridX;
    topology->gridY = gridY;
    int world_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    assert(world_size == gridX*gridY);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    topology->rankX = rank % gridX;
    topology->rankY = rank / gridX;

    if(topology->rankY==0) 
        MPI_Comm_split(MPI_COMM_WORLD, topology->rankY, topology->rankX, &(topology->scatterXComm));
    else
        MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, topology->rankX, &(topology->scatterXComm));
    MPI_Barrier(MPI_COMM_WORLD);

    if(topology->rankX==0)
        MPI_Comm_split(MPI_COMM_WORLD, topology->rankX, topology->rankY, &(topology->scatterYComm));
    else
        MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, topology->rankY, &(topology->scatterYComm));
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Comm_split(MPI_COMM_WORLD, topology->rankX, topology->rankY, &(topology->gatherXComm));
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Comm_split(MPI_COMM_WORLD, topology->rankY, topology->rankX, &(topology->gatherYComm));
    MPI_Barrier(MPI_COMM_WORLD);
}

void padding(particle_t** particles, int* nParticles, int minimalGrid) {
    int newSize = (*nParticles + minimalGrid - 1) / minimalGrid;
    newSize = newSize * minimalGrid;
    particle_t* oriParticles = *particles;
    particle_t* paddedParticles;
    paddedParticles = static_cast<particle_t*>(realloc(particles, newSize*sizeof(particle_t)));
    assert(paddedParticles != nullptr);
    if(paddedParticles != oriParticles) {
        particle_t* nppt = oriParticles + (*nParticles);
        for(int n=(*nParticles); n<newSize; n++) {
            nppt->enabled = false;
            nppt->ndim = oriParticles->ndim;
            nppt->nfeat = oriParticles->nfeat;
            nppt->id = UINT64_MAX-n;
            nppt->updateAcceleration = oriParticles->updateAcceleration;
            nppt++;
        }
    }
    *nParticles = newSize;
    *particles = paddedParticles;

    double* pos = oriParticles->position;
    double* vel = oriParticles->velocity;
    double* acc = oriParticles->acceleration;
    double* feat = oriParticles->features;

    double* paddedPos = static_cast<double*>(realloc(oriParticles->position, newSize*oriParticles->ndim*sizeof(double)));
    double* paddedVel = static_cast<double*>(realloc(oriParticles->velocity, newSize*oriParticles->ndim*sizeof(double)));
    double* paddedAcc = static_cast<double*>(realloc(oriParticles->acceleration, newSize*oriParticles->ndim*sizeof(double)));
    double* paddedFeatures = static_cast<double*>(realloc(oriParticles->features, newSize*oriParticles->nfeat*sizeof(double)));
    assert(paddedPos != nullptr);
    assert(paddedVel != nullptr);
    assert(paddedAcc != nullptr);
    assert(paddedFeatures != nullptr);

    particle_t* ppt = nullptr;
    double* dpt = nullptr;
    if(pos != paddedPos) {
        ppt = paddedParticles;
        dpt = paddedPos;
        for(int n=0; n<newSize; n++) {
            ppt->position = dpt;
            dpt += ppt->ndim;
            ppt++;
        }
    }
    if(vel != paddedVel) {
        ppt = paddedParticles;
        dpt = paddedVel;
        for(int n=0; n<newSize; n++) {
            ppt->velocity = dpt;
            dpt += ppt->ndim;
            ppt++;
        }
    }
    if(acc != paddedAcc) {
        ppt = paddedParticles;
        dpt = paddedAcc;
        for(int n=0; n<newSize; n++) {
            ppt->position = dpt;
            dpt += ppt->ndim;
            ppt++;
        }
    }
    if(feat != paddedFeatures) {
        ppt = paddedParticles;
        dpt = paddedFeatures;
        for(int n=0; n<newSize; n++) {
            ppt->features = dpt;
            dpt += ppt->nfeat;
            ppt++;
        }
    }
    
}

void padding(chunk_particles_t* chunkParticles, int minimalGrid) {
    particle_t* particles = chunkParticles->particles;
    int nParticles = chunkParticles->nParticle;
    padding(&particles, &nParticles, minimalGrid);
    chunkParticles->particles = particles;
    chunkParticles->nParticle = nParticles;
}

void shrink(particle_t** particles, int* nParticles) {
    particle_t* oriParticles = *particles;
    int oriNParticles = *nParticles;
    int newSize;

    particle_t* ppt = oriParticles;
    for(newSize=0; newSize<oriNParticles; newSize++) {
        if (!ppt->enabled) {
            break;
        }
        ppt++;
    }
    double* oriDpt; 
    double* newDpt; 

    oriDpt = oriParticles->position;
    newDpt = static_cast<double*>(realloc(oriDpt, newSize*oriParticles->ndim*sizeof(double)));
    assert(newDpt == oriDpt);

    oriDpt = oriParticles->velocity;
    newDpt = static_cast<double*>(realloc(oriDpt, newSize*oriParticles->ndim*sizeof(double)));
    assert(newDpt == oriDpt);
    
    oriDpt = oriParticles->acceleration;
    newDpt = static_cast<double*>(realloc(oriDpt, newSize*oriParticles->ndim*sizeof(double)));
    assert(newDpt == oriDpt);

    oriDpt = oriParticles->features;
    newDpt = static_cast<double*>(realloc(oriDpt, newSize*oriParticles->nfeat*sizeof(double)));
    assert(newDpt == oriDpt);

    ppt = static_cast<particle_t*>(realloc(oriParticles, newSize*sizeof(particle_t)));
    assert(ppt == oriParticles);

    *particles = ppt;
    *nParticles = newSize;
}

void shrink(chunk_particles_t* chunkParticles) {
    particle_t* particles = chunkParticles->particles;
    int nParticle = chunkParticles->nParticle;
    shrink(&particles, &nParticle);
    chunkParticles->particles = particles;
    chunkParticles->nParticle = nParticle;
}

