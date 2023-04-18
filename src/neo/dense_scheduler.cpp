#include "neo/dense_scheduler.h"
#include "neo/particle.h"
#include <vector>
#include <cassert>
#include <algorithm>

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

int minimalPaddingSize(int oriSize, int gridX, int gridY, int kernelSize) {
    int gridSize = gridX*gridY*kernelSize;
    int paddedSize = (oriSize+gridSize-1) / gridSize * gridSize;
    return paddedSize;
}

void padding(std::vector<particle_t>* particles, int paddedSize) {
    particles->resize(paddedSize);
    for(int n=particles->size(); n<paddedSize; n++) {
        particle_t& particle = particles->at(n);
        particle.enabled = false;
    }
}

void shrink(std::vector<particle_t>* particles) {
    particles->erase(
        std::remove_if(
            particles->begin(), 
            particles->end(), 
            [&](const particle_t& particle) {
                return !particle.enabled;
            }
        ), particles->end()
    );
}

void scatter(
    std::vector<particle_t>* localParticlesX, 
    std::vector<particle_t>* localParticlesY,
    std::vector<particle_t>* globalParticles,
    topology_t* topology) {
    const int& rankX = topology->rankX;
    const int& rankY = topology->rankY;
    const int& gridX = topology->gridX;
    const int& gridY = topology->gridY;
    if(rankX==0 && rankY==0) {
        assert(globalParticles->size()%gridX == 0);
        assert(globalParticles->size()%gridY == 0);
    }
    localParticlesX->clear();
    localParticlesY->clear();

    int* metadata = static_cast<int*>(malloc(2*sizeof(int)));
    assert(metadata != nullptr);
    if(rankX==0 && rankY==0) {
        *(metadata+0) = globalParticles->size() / gridX;
        *(metadata+1) = globalParticles->size() / gridY;
    }
    MPI_Bcast(metadata, 2, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    int localNParticlesX = *(metadata+0);
    int localNParticlesY = *(metadata+1);
    localParticlesX->resize(localNParticlesX);
    localParticlesY->resize(localNParticlesY);

    if(rankY==0) {
        if(rankX==0) {
            MPI_Scatter(
                const_cast<particle_t*>(globalParticles->data()),
                globalParticles->size() * sizeof(particle_t) / gridX,
                MPI_BYTE,
                const_cast<particle_t*>(localParticlesX->data()),
                localParticlesX->size() * sizeof(particle_t),
                MPI_BYTE,
                0,
                topology->scatterXComm
            );
        }
        else {
            MPI_Scatter(
                nullptr,
                -1,
                MPI_BYTE,
                const_cast<particle_t*>(localParticlesX->data()),
                localParticlesX->size() * sizeof(particle_t),
                MPI_BYTE,
                0,
                topology->scatterXComm
            );
        }
    }
    if(rankX==0) {
        if(rankY==0) {
            MPI_Scatter(
                const_cast<particle_t*>(globalParticles->data()),
                globalParticles->size() * sizeof(particle_t) / gridY,
                MPI_BYTE,
                const_cast<particle_t*>(localParticlesY->data()),
                localParticlesY->size() * sizeof(particle_t),
                MPI_BYTE,
                0,
                topology->scatterYComm
            );
        }
        else {
            MPI_Scatter(
                nullptr,
                -1,
                MPI_BYTE,
                const_cast<particle_t*>(localParticlesY->data()),
                localParticlesY->size() * sizeof(particle_t),
                MPI_BYTE,
                0,
                topology->scatterYComm
            );
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(
        const_cast<particle_t*>(localParticlesX->data()),
        localParticlesX->size() * sizeof(particle_t),
        MPI_BYTE,
        0,
        topology->reduceXComm
    );
    MPI_Bcast(
        const_cast<particle_t*>(localParticlesY->data()),
        localParticlesY->size() * sizeof(particle_t),
        MPI_BYTE,
        0,
        topology->reduceYComm
    );
    MPI_Barrier(MPI_COMM_WORLD);
}

void reduce(
    std::vector<particle_t>* localParticlesX,
    std::vector<particle_t>* localParticlesY,
    topology_t* topology) {
    double *accelerationsXTx, *accelerationsXRx;
    double *accelerationsYTx, *accelerationsYRx;
    accelerationsXTx = static_cast<double*>(
        malloc(localParticlesX->size() * sizeof(double) * PARTICLE_N_DIM));
    accelerationsYTx = static_cast<double*>(
        malloc(localParticlesY->size() * sizeof(double) * PARTICLE_N_DIM));
    accelerationsXRx = static_cast<double*>(
        malloc(localParticlesX->size() * sizeof(double) * PARTICLE_N_DIM));
    accelerationsYRx = static_cast<double*>(
        malloc(localParticlesY->size() * sizeof(double) * PARTICLE_N_DIM));
    
    #pragma omp parallel for collapse(2)
    for(int n=0; n<localParticlesX->size(); n++) {
        for(int nn=0; nn<PARTICLE_N_DIM; nn++) {
            *(accelerationsXTx+n*PARTICLE_N_DIM+nn) = 
                localParticlesX->at(n).acceleration[nn];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(
        accelerationsXTx,
        accelerationsXRx,
        localParticlesX->size() * PARTICLE_N_DIM,
        MPI_DOUBLE,
        MPI_SUM,
        topology->reduceXComm
    );
    #pragma omp parallel for collapse(2)
    for(int n=0; n<localParticlesX->size(); n++) {
        for(int nn=0; nn<PARTICLE_N_DIM; nn++) {
            localParticlesX->at(n).acceleration[nn] =
                *(accelerationsXRx+n*PARTICLE_N_DIM+nn);
        }
    }


    #pragma omp parallel for collapse(2)
    for(int n=0; n<localParticlesY->size(); n++) {
        for(int nn=0; nn<PARTICLE_N_DIM; nn++) {
            *(accelerationsYTx+n*PARTICLE_N_DIM+nn) = 
                localParticlesY->at(n).acceleration[nn];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(
        accelerationsYTx,
        accelerationsYRx,
        localParticlesY->size() * PARTICLE_N_DIM,
        MPI_DOUBLE,
        MPI_SUM,
        topology->reduceYComm
    );
    #pragma omp parallel for collapse(2)
    for(int n=0; n<localParticlesY->size(); n++) {
        for(int nn=0; nn<PARTICLE_N_DIM; nn++) {
            localParticlesY->at(n).acceleration[nn] =
                *(accelerationsYRx+n*PARTICLE_N_DIM+nn);
        }
    }
    free(accelerationsXTx);
    free(accelerationsXRx);
    free(accelerationsYTx);
    free(accelerationsYRx);
}

void gather(
    std::vector<particle_t>* localParticlesX,
    std::vector<particle_t>* localParticlesY,
    std::vector<particle_t>* globalParticles,
    topology_t* topology) {
    const int& rankX = topology->rankX;
    const int& rankY = topology->rankY;
    const int& gridX = topology->gridX;
    const int& gridY = topology->gridY;

    if(rankX==0 && rankY==0) {
        globalParticles->resize(localParticlesX->size()*gridX);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(rankY==0) {
        if(rankX==0) {
            MPI_Gather(
                const_cast<particle_t*>(localParticlesX->data()),
                localParticlesX->size() * sizeof(particle_t),
                MPI_BYTE,
                const_cast<particle_t*>(globalParticles->data()),
                globalParticles->size() * sizeof(particle_t) / gridX,
                MPI_BYTE,
                0,
                topology->scatterXComm
            );
        } 
        else {
            MPI_Gather(
                const_cast<particle_t*>(localParticlesX->data()),
                localParticlesX->size() * sizeof(particle_t),
                MPI_BYTE,
                nullptr,
                -1,
                MPI_BYTE,
                0,
                topology->scatterXComm
            );
        }
    }
}

void update(
    std::vector<particle_t>* localParticlesX,
    std::vector<particle_t>* localParticlesY,
    topology_t* topology,
    double timeStep,
    void (*updateAccelerationFn)(particle_t*, particle_t*)) {

    // clear Accelerations
    for (int n=0; n<localParticlesX->size(); n++) 
        clearAcceleration(std::addressof(localParticlesX->at(n)));
    for (int n=0; n<localParticlesY->size(); n++) 
        clearAcceleration(std::addressof(localParticlesY->at(n)));
    
    // update Accelerations locally
    for (int nx=0; nx<localParticlesX->size(); nx++)
        for (int ny=0; ny<localParticlesY->size(); ny++) {
            particle_t* particleX = std::addressof(localParticlesX->at(nx));
            particle_t* particleY = std::addressof(localParticlesY->at(ny));
            if (particleX->id == particleY->id)
                continue;
            updateAccelerationFn(particleX, particleY);
        }
    
    // update Accelerations globally
    reduce(localParticlesX, localParticlesY, topology);

    // update vel and pos
    for (int n=0; n<localParticlesX->size(); n++) {
        updateVelocity(std::addressof(localParticlesX->at(n)), timeStep);
        updatePosition(std::addressof(localParticlesX->at(n)), timeStep);
    }
    for (int n=0; n<localParticlesY->size(); n++) {
        updateVelocity(std::addressof(localParticlesY->at(n)), timeStep);
        updatePosition(std::addressof(localParticlesY->at(n)), timeStep);
    }
}