#include "dense_scheduler.h"
#include <mpi.h>
#include <stdlib.h>
#include <cassert>

__inline__ int sizeAfterPadding(int nParticles, int minimalGrid) {
    return (nParticles+minimalGrid-1)/minimalGrid * minimalGrid;
}

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

    MPI_Comm_split(MPI_COMM_WORLD, topology->rankX, topology->rankY, &(topology->reduceXComm));
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Comm_split(MPI_COMM_WORLD, topology->rankY, topology->rankX, &(topology->reduceYComm));
    MPI_Barrier(MPI_COMM_WORLD);
}

void padding(particle_t** particles, int* nParticles, int minimalGrid) {
    int oriSize = *nParticles;
    int newSize = sizeAfterPadding(oriSize, minimalGrid);
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

void scatter(particle_t** localParticlesX, int* localNParticlesX,
            particle_t** localParticlesY, int* localNParticlesY,
            particle_t* wholeParticles, int wholeNParticles, 
            topology_t* topology) {
    bool wholeParticlesInitNull = false;
    if(wholeParticles==nullptr) {
        wholeParticlesInitNull = true;
        alloc_particles(&wholeParticles, 1, 1, 1);
    }

    const int& rankX = topology->rankX;
    const int& rankY = topology->rankY;
    const int& gridX = topology->gridX;
    const int& gridY = topology->gridY;
    if(rankX==0 && rankY==0) {
        assert(wholeNParticles%gridX == 0);
        assert(wholeNParticles%gridY == 0);
    }

    int* localParticlesMetadata = static_cast<int*>(malloc(4*sizeof(int)));
    assert(localParticlesMetadata != nullptr);

    if(rankX==0 && rankY==0) {
        *(localParticlesMetadata+0) = wholeNParticles / gridX;
        *(localParticlesMetadata+1) = wholeNParticles / gridY;
        *(localParticlesMetadata+2) = wholeParticles->ndim;
        *(localParticlesMetadata+3) = wholeParticles->nfeat;
    }

    MPI_Bcast(localParticlesMetadata, 4, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    int ndim, nfeat;
    *localNParticlesX = *(localParticlesMetadata+0);
    *localNParticlesY = *(localParticlesMetadata+1);
    ndim = *(localParticlesMetadata+2);
    nfeat = *(localParticlesMetadata+3);
    alloc_particles(localParticlesX, *localNParticlesX, ndim, nfeat);
    alloc_particles(localParticlesY, *localNParticlesY, ndim, nfeat);

    double *pos, *vel, *acc, *feat;

    if(rankY==0) {
        pos = (*localParticlesX)->position;
        vel = (*localParticlesX)->velocity;
        acc = (*localParticlesX)->acceleration;
        feat = (*localParticlesX)->features;
        printf("%s:%d, %lx %lx %lx %lx\n", __FILE__, __LINE__, pos, vel, acc, feat);

        MPI_Barrier(topology->scatterXComm);
        printf("(%d, %d) %s %d %d@%d\n", topology->rankX, topology->rankY, __FILE__, __LINE__, *localNParticlesX, wholeNParticles);
        MPI_Scatter(
            wholeParticles,
            wholeNParticles*sizeof(particle_t)/gridX,
            MPI_BYTE,
            (*localParticlesX),
            (*localNParticlesX)*sizeof(particle_t),
            MPI_BYTE,
            0,
            topology->scatterXComm
        );

        #pragma omp parallel for
        for (int n=0; n<*localNParticlesX; n++) {
            ((*localParticlesX)+n)->position = pos+n*ndim;
            ((*localParticlesX)+n)->velocity = vel+n*ndim;
            ((*localParticlesX)+n)->acceleration = acc+n*ndim;
            ((*localParticlesX)+n)->features = feat+n*nfeat;
        }

        MPI_Barrier(topology->scatterXComm);
        MPI_Scatter(
            wholeParticles->position,
            wholeNParticles*ndim/gridX,
            MPI_DOUBLE,
            pos, 
            (*localNParticlesX)*ndim,
            MPI_DOUBLE,
            0,
            topology->scatterXComm
        );

        MPI_Barrier(topology->scatterXComm);
        printf("(%d, %d) %s %d %d@%d\n", topology->rankX, topology->rankY, __FILE__, __LINE__, *localNParticlesX, wholeNParticles);
        MPI_Scatter(
            wholeParticles->velocity,
            wholeNParticles*ndim/gridX,
            MPI_DOUBLE,
            (*localParticlesX)->velocity, 
            (*localNParticlesX)*ndim,
            MPI_DOUBLE,
            0,
            topology->scatterXComm
        );
        
        MPI_Barrier(topology->scatterXComm);
        printf("(%d, %d) %s %d %d@%d\n", topology->rankX, topology->rankY, __FILE__, __LINE__, *localNParticlesX, wholeNParticles);
        MPI_Scatter(
            wholeParticles->acceleration,
            wholeNParticles*ndim/gridX,
            MPI_DOUBLE,
            (*localParticlesX)->acceleration, 
            (*localNParticlesX)*ndim,
            MPI_DOUBLE,
            0,
            topology->scatterXComm
        );
        
        MPI_Barrier(topology->scatterXComm);
        printf("(%d, %d) %s %d %d@%d\n", topology->rankX, topology->rankY, __FILE__, __LINE__, *localNParticlesX, wholeNParticles);
        MPI_Scatter(
            wholeParticles->features,
            wholeNParticles*nfeat/gridX,
            MPI_DOUBLE,
            (*localParticlesX)->features, 
            (*localNParticlesX)*nfeat,
            MPI_DOUBLE,
            0,
            topology->scatterXComm
        );
    }

    pos = (*localParticlesX)->position;
    vel = (*localParticlesX)->velocity;
    acc = (*localParticlesX)->acceleration;
    feat = (*localParticlesX)->features;

    MPI_Bcast(
        (*localParticlesX),
        (*localNParticlesX)*sizeof(particle_t),
        MPI_BYTE,
        0,
        topology->reduceXComm
    );

    #pragma omp parallel for
    for (int n=0; n<*localNParticlesX; n++) {
        (*localParticlesX+n)->position = pos+n*ndim;
        (*localParticlesX+n)->velocity = vel+n*ndim;
        (*localParticlesX+n)->acceleration = acc+n*ndim;
        (*localParticlesX+n)->features = feat+n*nfeat;
    }
    MPI_Bcast(
        (*localParticlesX)->position,
        (*localNParticlesX)*ndim,
        MPI_DOUBLE,
        0,
        topology->reduceXComm
    );
    MPI_Bcast(
        (*localParticlesX)->velocity,
        (*localNParticlesX)*ndim,
        MPI_DOUBLE,
        0,
        topology->reduceXComm
    );
    MPI_Bcast(
        (*localParticlesX)->acceleration,
        (*localNParticlesX)*ndim,
        MPI_DOUBLE,
        0,
        topology->reduceXComm
    );
    MPI_Bcast(
        (*localParticlesX)->features,
        (*localNParticlesX)*nfeat,
        MPI_DOUBLE,
        0,
        topology->reduceXComm
    );


    if(rankX==0) {
        pos = (*localParticlesY)->position;
        vel = (*localParticlesY)->velocity;
        acc = (*localParticlesY)->acceleration;
        feat = (*localParticlesY)->features;

        MPI_Scatter(
            wholeParticles,
            wholeNParticles*sizeof(particle_t)/gridY,
            MPI_BYTE,
            (*localParticlesY),
            (*localNParticlesY)*sizeof(particle_t),
            MPI_BYTE,
            0,
            topology->scatterYComm
        );

        #pragma omp parallel for
        for (int n=0; n<*localNParticlesY; n++) {
            (*localParticlesY+n)->position = pos+n*ndim;
            (*localParticlesY+n)->velocity = vel+n*ndim;
            (*localParticlesY+n)->acceleration = acc+n*ndim;
            (*localParticlesY+n)->features = feat+n*nfeat;
        }

        MPI_Scatter(
            wholeParticles->position,
            wholeNParticles*ndim/gridY,
            MPI_DOUBLE,
            (*localParticlesY)->position, 
            (*localNParticlesY)*ndim,
            MPI_DOUBLE,
            0,
            topology->scatterYComm
        );
        MPI_Scatter(
            wholeParticles->velocity,
            wholeNParticles*ndim/gridY,
            MPI_DOUBLE,
            (*localParticlesY)->velocity, 
            (*localNParticlesY)*ndim,
            MPI_DOUBLE,
            0,
            topology->scatterYComm
        );
        MPI_Scatter(
            wholeParticles->acceleration,
            wholeNParticles*ndim/gridY,
            MPI_DOUBLE,
            (*localParticlesY)->acceleration, 
            (*localNParticlesY)*ndim,
            MPI_DOUBLE,
            0,
            topology->scatterYComm
        );
        MPI_Scatter(
            wholeParticles->features,
            wholeNParticles*nfeat/gridY,
            MPI_DOUBLE,
            (*localParticlesY)->features, 
            (*localNParticlesY)*nfeat,
            MPI_DOUBLE,
            0,
            topology->scatterYComm
        );
    }

    pos = (*localParticlesY)->position;
    vel = (*localParticlesY)->velocity;
    acc = (*localParticlesY)->acceleration;
    feat = (*localParticlesY)->features;

    MPI_Bcast(
        (*localParticlesY),
        (*localNParticlesY)*sizeof(particle_t),
        MPI_BYTE,
        0,
        topology->reduceYComm
    );

    #pragma omp parallel for
    for (int n=0; n<*localNParticlesY; n++) {
        (*localParticlesY+n)->position = pos+n*ndim;
        (*localParticlesY+n)->velocity = vel+n*ndim;
        (*localParticlesY+n)->acceleration = acc+n*ndim;
        (*localParticlesY+n)->features = feat+n*nfeat;
    }
    MPI_Bcast(
        (*localParticlesY)->position,
        (*localNParticlesY)*ndim,
        MPI_DOUBLE,
        0,
        topology->reduceYComm
    );
    MPI_Bcast(
        (*localParticlesY)->velocity,
        (*localNParticlesY)*ndim,
        MPI_DOUBLE,
        0,
        topology->reduceYComm
    );
    MPI_Bcast(
        (*localParticlesY)->acceleration,
        (*localNParticlesY)*ndim,
        MPI_DOUBLE,
        0,
        topology->reduceYComm
    );
    MPI_Bcast(
        (*localParticlesY)->features,
        (*localNParticlesY)*nfeat,
        MPI_DOUBLE,
        0,
        topology->reduceYComm
    );
    if(wholeParticlesInitNull) {
        free_particles(wholeParticles);
        wholeParticles = nullptr;
    }
}

void scatter(chunk_particles_t** localChunkParticlesX, 
                chunk_particles_t** localChunkParticlesY,
                chunk_particles_t* wholeChunkParticles,
                topology_t* topology) {
    *localChunkParticlesX = static_cast<chunk_particles_t*>(malloc(sizeof(chunk_particles_t)));
    *localChunkParticlesY = static_cast<chunk_particles_t*>(malloc(sizeof(chunk_particles_t)));
    particle_t *localParticlesX=nullptr, *localParticlesY=nullptr, *wholeParticles=nullptr;
    int localNParticlesX=-1, localNParticlesY=-1, wholeNParticles;
    if(wholeChunkParticles==nullptr){
        wholeParticles = nullptr;
        wholeNParticles = -1;
    }
    else {
        wholeParticles = wholeChunkParticles->particles;
        wholeNParticles = wholeChunkParticles->nParticle;
    }
    scatter(
        &localParticlesX, &localNParticlesX, &localParticlesY, &localNParticlesY,
        wholeParticles, wholeNParticles, topology
    );
    (*localChunkParticlesX)->particles = localParticlesX;
    (*localChunkParticlesX)->nParticle = localNParticlesX;
    (*localChunkParticlesY)->particles = localParticlesY;
    (*localChunkParticlesY)->nParticle = localNParticlesY;
}

void reduce(particle_t* localParticlesX, int localNParticlesX, 
        particle_t* localParticlesY, int localNParticlesY, topology_t* topology) {
    double* recvBuf1;
    recvBuf1 = static_cast<double*>(malloc(localNParticlesX*localParticlesX->ndim*sizeof(double)));
    MPI_Allreduce(
        localParticlesX->acceleration,
        recvBuf1,
        localNParticlesX*localParticlesX->ndim,
        MPI_DOUBLE,
        MPI_SUM,
        topology->reduceXComm
    );
    #pragma omp parallel for collapse(2)
    for(int n=0; n<localNParticlesX; n++)
        for(int nn=0; nn<localParticlesX->ndim; nn++)
            *((localParticlesX->acceleration) + n*localParticlesX->ndim + nn) = 
                    *(recvBuf1 + n*localParticlesX->ndim + nn);
    
    double* recvBuf2;
    recvBuf2 = static_cast<double*>(malloc(localNParticlesY*(localParticlesY->ndim)*sizeof(double)));
    MPI_Allreduce(
        localParticlesY->acceleration,
        recvBuf2,
        localNParticlesY*localParticlesY->ndim,
        MPI_DOUBLE,
        MPI_SUM,
        topology->reduceYComm
    );
    #pragma omp parallel for collapse(2)
    for(int n=0; n<localNParticlesY; n++)
        for(int nn=0; nn<localParticlesY->ndim; nn++)
            *((localParticlesY->acceleration)+n*localParticlesY->ndim+nn) = 
                    *(recvBuf2+n*localParticlesY->ndim+nn);
    free(recvBuf1);
    free(recvBuf2);
}

void reduce(chunk_particles_t* localChunkParticlesX, chunk_particles_t* localChunkParticlesY, topology_t* topology) {
    reduce(localChunkParticlesX->particles, localChunkParticlesX->nParticle, 
            localChunkParticlesY->particles, localChunkParticlesY->nParticle, topology);
}


void gather(particle_t* localParticlesX, int localNParticlesX, 
                particle_t* localParticlesY, int localNParticlesY,
                particle_t** globalParticles, int* globalNParticles, topology_t* topology) {
    const int& rankX = topology->rankX;
    const int& rankY = topology->rankY;
    const int& gridX = topology->gridX;
    const int& gridY = topology->gridY;

    if(rankX==0 && rankY==0) {
        assert(*globalNParticles == localNParticlesX * gridX);
        assert(*globalNParticles == localNParticlesY * gridY);
    }

    double *pos, *vel, *acc, *feat;
    if(rankX==0) {
        if(rankY==0) {
            pos = (*globalParticles)->position;
            vel = (*globalParticles)->velocity;
            acc = (*globalParticles)->acceleration;
            feat = (*globalParticles)->features;

            MPI_Gather(localParticlesY, 
                localNParticlesY*sizeof(particle_t),
                MPI_BYTE,
                *globalParticles,
                *globalNParticles*sizeof(particle_t)/gridY,
                MPI_BYTE,
                0, 
                topology->scatterYComm
            );

            MPI_Gather(localParticlesY->position,
                localNParticlesY*localParticlesY->ndim,
                MPI_DOUBLE,
                pos,
                (*globalNParticles)*(*globalParticles)->ndim/gridY,
                MPI_DOUBLE,
                0,
                topology->scatterYComm
            );

            MPI_Gather(localParticlesY->velocity,
                localNParticlesY*localParticlesY->ndim,
                MPI_DOUBLE,
                vel,
                (*globalNParticles)*(*globalParticles)->ndim/gridY,
                MPI_DOUBLE,
                0,
                topology->scatterYComm
            );

            MPI_Gather(localParticlesY->acceleration,
                localNParticlesY*localParticlesY->ndim,
                MPI_DOUBLE,
                acc,
                (*globalNParticles)*(*globalParticles)->ndim/gridY,
                MPI_DOUBLE,
                0,
                topology->scatterYComm
            );

            MPI_Gather(localParticlesY->features,
                localNParticlesY*localParticlesY->nfeat,
                MPI_DOUBLE,
                feat,
                (*globalNParticles)*(*globalParticles)->nfeat/gridY,
                MPI_DOUBLE,
                0,
                topology->scatterYComm
            );
        }
        else {
            MPI_Gather(localParticlesY, 
                localNParticlesY*sizeof(particle_t),
                MPI_BYTE,
                nullptr,
                -1,
                MPI_BYTE,
                0, 
                topology->scatterYComm
            );

            MPI_Gather(localParticlesY->position,
                localNParticlesY*localParticlesY->ndim,
                MPI_DOUBLE,
                nullptr,
                -1,
                MPI_DOUBLE,
                0,
                topology->scatterYComm
            );

            MPI_Gather(localParticlesY->velocity,
                localNParticlesY*localParticlesY->ndim,
                MPI_DOUBLE,
                nullptr,
                -1,
                MPI_DOUBLE,
                0,
                topology->scatterYComm
            );

            MPI_Gather(localParticlesY->acceleration,
                localNParticlesY*localParticlesY->ndim,
                MPI_DOUBLE,
                nullptr,
                -1,
                MPI_DOUBLE,
                0,
                topology->scatterYComm
            );

            MPI_Gather(localParticlesY->features,
                localNParticlesY*localParticlesY->nfeat,
                MPI_DOUBLE,
                nullptr,
                -1,
                MPI_DOUBLE,
                0,
                topology->scatterYComm
            );
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(rankX==0 && rankY==0) {
        #pragma omp parallel for
        for(int n=0; n<(*globalNParticles); n++) {
            ((*globalParticles)+n)->position = pos+n*((*globalParticles)->ndim);
            ((*globalParticles)+n)->velocity = vel+n*((*globalParticles)->ndim);
            ((*globalParticles)+n)->acceleration = acc+n*((*globalParticles)->ndim);
            ((*globalParticles)+n)->features = feat+n*((*globalParticles)->nfeat);
        }
    }
    return;
}

void gather(chunk_particles_t* localChunkParticlesX, 
            chunk_particles_t* localChunkParticlesY, 
            chunk_particles_t* globalChunkParticles, 
            topology_t* topology) {
    particle_t* localParticlesX = localChunkParticlesX->particles;
    int localNParticlesX = localChunkParticlesX->nParticle;
    particle_t* localParticlesY = localChunkParticlesY->particles;
    int localNParticlesY = localChunkParticlesY->nParticle;
    particle_t* globalParticles;
    int globalNParticles;
    if(globalChunkParticles==nullptr) {
        globalParticles = nullptr;
        globalNParticles = -1;
    }
    else {
        globalParticles = globalChunkParticles->particles;
        globalNParticles = globalChunkParticles->nParticle;
        printf("%s:%d %lx %d\n", __FILE__, __LINE__, globalParticles, globalNParticles);
    }

    gather(localParticlesX, localNParticlesX, localParticlesY, 
            localNParticlesY, &globalParticles, &globalNParticles, topology);
    if(globalChunkParticles!=nullptr) {
        printf("%s:%d %lx %d\n", __FILE__, __LINE__, globalParticles, globalNParticles);
        globalChunkParticles->particles = globalParticles;
        globalChunkParticles->nParticle = globalNParticles;
    }
}

void update(particle_t* particlesX, int nParticlesX, particle_t* particlesY, int nParticlesY, topology_t* topology, double timeStep) {
    // clear accelerations
    #pragma omp parallel for collapse(2)
    for(int nx=0; nx<nParticlesX; nx++)
        for(int nn=0; nn<particlesX->ndim; nn++)
            *((particlesX+nx)->acceleration + nn) = 0.0;
    #pragma omp parallel for collapse(2)
    for(int ny=0; ny<nParticlesY; ny++)
        for(int nn=0; nn<particlesY->ndim; nn++)
            *((particlesY+ny)->acceleration + nn) = 0.0;
    
    // update accelerations local
    #pragma omp parallel for collapse(2)
    for(int nx=0; nx<nParticlesX; nx++)
        for(int ny=0; ny<nParticlesY; ny++) 
            if((particlesX+nx)->id == (particlesY+nx)->id)
                updateAcceleration(particlesX+nx, particlesY+ny);
    
    // update accelerations global
    reduce(particlesX, nParticlesX, particlesY, nParticlesY, topology);

    // update velocity and position
    #pragma omp parallel for collapse(2)
    for(int nx=0; nx<nParticlesX; nx++)
        for(int ndim=0; ndim<particlesX->ndim; ndim++) {
            *(particlesX->velocity+ndim) += timeStep * *(particlesX->acceleration+ndim);
            *(particlesX->position+ndim) += timeStep * *(particlesX->velocity+ndim);
        }
    #pragma omp parallel for collapse(2)
    for(int ny=0; ny<nParticlesY; ny++)
        for(int ndim=0; ndim<particlesY->ndim; ndim++) {
            *(particlesY->velocity+ndim) += timeStep * *(particlesY->acceleration+ndim);
            *(particlesY->position+ndim) += timeStep * *(particlesY->velocity+ndim);
        }
}

void update(chunk_particles_t* chunkParticlesX, 
            chunk_particles_t* chunkParticlesY, 
            topology_t* topology, double timeStep) {
    update(chunkParticlesX->particles, chunkParticlesX->nParticle, 
            chunkParticlesY->particles, chunkParticlesY->nParticle, topology, timeStep);
}