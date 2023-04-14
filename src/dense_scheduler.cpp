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
    int oriSize = *nParticles;
    int newSize = (*nParticles + minimalGrid - 1) / minimalGrid * minimalGrid;
    particle_t *oriP, *newP;
    double *oriPos, *oriVel, *oriAcc, *oriFeat;
    double *newPos, *newVel, *newAcc, *newFeat;
    oriP = *particles;
    oriPos = oriP->position;
    oriVel = oriP->velocity;
    oriAcc = oriP->acceleration;
    oriFeat = oriP->features;

    int ndim = oriP->ndim;
    int nfeat = oriP->nfeat;

    newP = static_cast<particle_t*>(realloc(oriP, newSize*sizeof(particle_t)));
    newPos = static_cast<double*>(realloc(oriPos, newSize*ndim*sizeof(double)));
    newVel = static_cast<double*>(realloc(oriVel, newSize*ndim*sizeof(double)));
    newAcc = static_cast<double*>(realloc(oriAcc, newSize*ndim*sizeof(double)));
    newFeat = static_cast<double*>(realloc(oriFeat, newSize*nfeat*sizeof(double)));

    assert(newP != nullptr);    assert(newPos != nullptr);      
    assert(newVel != nullptr);  assert(newAcc != nullptr);

    #pragma omp parallel for collapse(2)
    for(int n=oriSize; n<newSize; n++) 
        for(int nn=0; nn<ndim; nn++) {
            *(newPos+n*ndim+nn) = 0;
            *(newVel+n*ndim+nn) = 0;
            *(newAcc+n*ndim+nn) = 0;
        }
    
    #pragma omp parallel for collapse(2)
    for(int n=oriSize; n<newSize; n++)
        for(int nn=0; nn<nfeat; nn++)
            *(newFeat+n*nfeat+nn) = 0;
    
    #pragma omp parallel for
    for(int n=0; n<newSize; n++) {
        if((newPos != oriPos) || (n>=oriSize))
            (newP+n)->position = newPos+n*ndim;
        if((newVel != oriVel) || (n>=oriSize))
            (newP+n)->velocity = newVel+n*ndim;
        if((newAcc != oriAcc) || (n>=oriSize))
            (newP+n)->acceleration = newAcc+n*ndim;
        if((newFeat != oriFeat) || (n>=oriSize))
            (newP+n)->features = newFeat+n*nfeat;
        if(n>=oriSize) {
            (newP+n)->enabled = false;
            (newP+n)->id = UINT64_MAX - n;
            (newP+n)->ndim = ndim;
            (newP+n)->nfeat = nfeat;
            (newP+n)->updateAcceleration = newP->updateAcceleration;
        }
    }

    *particles = newP;
    *nParticles = newSize;
}


// void padding(particle_t** particles, int* nParticles, int minimalGrid) {
//     int newSize = (*nParticles + minimalGrid - 1) / minimalGrid;
//     newSize = newSize * minimalGrid;
//     particle_t* oriParticles = *particles;
//     particle_t* paddedParticles;
//     paddedParticles = static_cast<particle_t*>(realloc(*particles, newSize*sizeof(particle_t)));
//     assert(paddedParticles != nullptr);
//     if(paddedParticles != oriParticles) {
//         particle_t* nppt = paddedParticles + (*nParticles);
//         for(int n=(*nParticles); n<newSize; n++) {
//             nppt->enabled = false;
//             nppt->ndim = paddedParticles->ndim;
//             nppt->nfeat = paddedParticles->nfeat;
//             nppt->id = UINT64_MAX-n;
//             nppt->updateAcceleration = paddedParticles->updateAcceleration;
//             nppt++;
//         }
//     }
//     *nParticles = newSize;
//     *particles = paddedParticles;

//     double* pos = paddedParticles->position;
//     double* vel = paddedParticles->velocity;
//     double* acc = paddedParticles->acceleration;
//     double* feat = paddedParticles->features;

//     double* paddedPos = static_cast<double*>(realloc(paddedParticles->position, newSize*paddedParticles->ndim*sizeof(double)));
//     double* paddedVel = static_cast<double*>(realloc(paddedParticles->velocity, newSize*paddedParticles->ndim*sizeof(double)));
//     double* paddedAcc = static_cast<double*>(realloc(paddedParticles->acceleration, newSize*paddedParticles->ndim*sizeof(double)));
//     double* paddedFeatures = static_cast<double*>(realloc(paddedParticles->features, newSize*paddedParticles->nfeat*sizeof(double)));
//     assert(paddedPos != nullptr);
//     assert(paddedVel != nullptr);
//     assert(paddedAcc != nullptr);
//     assert(paddedFeatures != nullptr);

//     particle_t* ppt = nullptr;
//     double* dpt = nullptr;
//     if(pos != paddedPos) {
//         ppt = paddedParticles;
//         dpt = paddedPos;
//         for(int n=0; n<newSize; n++) {
//             ppt->position = dpt;
//             for(int nn=0; nn<ppt->ndim; nn++)
//                 *(dpt++) = 0.0;
//             // dpt += ppt->ndim;
//             ppt++;
//         }
//     }
//     if(vel != paddedVel) {
//         ppt = paddedParticles;
//         dpt = paddedVel;
//         for(int n=0; n<newSize; n++) {
//             ppt->velocity = dpt;
//             for(int nn=0; nn<ppt->ndim; nn++)
//                 *(dpt++) = 0.0;
//             // dpt += ppt->ndim;
//             ppt++;
//         }
//     }
//     if(acc != paddedAcc) {
//         ppt = paddedParticles;
//         dpt = paddedAcc;
//         for(int n=0; n<newSize; n++) {
//             ppt->acceleration = dpt;
//             for(int nn=0; nn<ppt->ndim; nn++)
//                 *(dpt++) = 0.0;
//                 // dpt += ppt->ndim;
//             ppt++;
//         }
//     }
//     if(feat != paddedFeatures) {
//         ppt = paddedParticles;
//         dpt = paddedFeatures;
//         for(int n=0; n<newSize; n++) {
//             ppt->features = dpt;
//             for(int nn=0; nn<ppt->nfeat; nn++)
//                 *(dpt++) = 0.0;
//             ppt++;
//         }
//     }
    
// }

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
    double* dpt;
    for(newSize=0; newSize<oriNParticles; newSize++) {
        if (!ppt->enabled) {
            break;
        }
        ppt++;
    }
    double *oriPos, *oriVel, *oriAcc, *oriFeat;
    double *newPos, *newVel, *newAcc, *newFeat; 
    particle_t *oriP, *newP;

    oriPos = oriParticles->position;
    newPos = static_cast<double*>(realloc(oriPos, newSize*oriParticles->ndim*sizeof(double)));
    assert(newPos != nullptr);

    oriVel = oriParticles->velocity;
    newVel = static_cast<double*>(realloc(oriVel, newSize*oriParticles->ndim*sizeof(double)));
    assert(newVel != nullptr);
    
    oriAcc = oriParticles->acceleration;
    newAcc = static_cast<double*>(realloc(oriAcc, newSize*oriParticles->ndim*sizeof(double)));
    assert(newAcc != nullptr);

    oriFeat = oriParticles->features;
    newFeat = static_cast<double*>(realloc(oriFeat, newSize*oriParticles->nfeat*sizeof(double)));
    assert(newFeat != nullptr);

    oriP = oriParticles;
    newP = static_cast<particle_t*>(realloc(oriP, newSize*sizeof(particle_t)));
    assert(newP != nullptr);

    if(oriPos != newPos) {
        ppt = newP;
        dpt = newPos;
        for(int n=0; n<newSize; n++) {
            ppt->position = dpt;
            dpt+=ppt->ndim;
            ppt++;
        }
    }

    if(oriVel != newVel) {
        ppt = newP;
        dpt = newVel;
        for(int n=0; n<newSize; n++) {
            ppt->velocity = dpt;
            dpt+=ppt->ndim;
            ppt++;
        }
    }

    if(oriAcc != newAcc) {
        ppt = newP;
        dpt = newAcc;
        for(int n=0; n<newSize; n++) {
            ppt->acceleration = dpt;
            dpt+=ppt->ndim;
            ppt++;
        }
    }

    if(oriFeat != newFeat) {
        ppt = newP;
        dpt = newFeat;
        for(int n=0; n<newSize; n++) {
            ppt->features = dpt;
            dpt+=ppt->nfeat;
            ppt++;
        }
    }

    *particles = newP;
    *nParticles = newSize;
}

void shrink(chunk_particles_t* chunkParticles) {
    particle_t* particles = chunkParticles->particles;
    int nParticle = chunkParticles->nParticle;
    shrink(&particles, &nParticle);
    chunkParticles->particles = particles;
    chunkParticles->nParticle = nParticle;
}

