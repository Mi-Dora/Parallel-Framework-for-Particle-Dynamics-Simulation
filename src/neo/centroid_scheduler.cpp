#include "neo/centroid_scheduler.h"
#include "neo/particle.h"
#include <vector>
#include <algorithm>
#include <tuple>
#include <mpi.h>
#include <cassert>

class _CompParticlesPos {
public:
    bool operator() (std::tuple<particle_t*, double> a, std::tuple<particle_t*, double> b) {
        return (a<b);
    }
};

void setupCommunicators(topology_t* topology, int gridX, int gridY, int gridZ, int subGridX, int subGridY, int subGridZ) {
    topology->gridX = gridX;
    topology->gridY = gridY;
    topology->gridZ = gridZ;
    topology->subgridX = subGridX;
    topology->subgridY = subGridY;
    topology->subgridZ = subGridZ;
    int compSize, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &compSize);
    topology->rankX = rank % gridX;
    topology->rankY = (rank / gridX) % gridY;
    topology->rankZ = rank / gridX / gridY;

    MPI_Comm_split(MPI_COMM_WORLD, topology->rankY+topology->rankZ*gridY, topology->rankX, &(topology->xEdgeComm));
    MPI_Comm_split(MPI_COMM_WORLD, topology->rankX+topology->rankZ*gridX, topology->rankY, &(topology->yEdgeComm));
    MPI_Comm_split(MPI_COMM_WORLD, topology->rankY+topology->rankX*gridY, topology->rankZ, &(topology->zEdgeComm));
}


void divideGroupN(std::vector<particle_t>* allParticles, 
                std::vector<std::vector<particle_t>>* groupParticles,
                int ndim,
                int ngroup) {
    std::vector<std::tuple<particle_t*, double>> idPosPairs;
    for(int n=0; n<allParticles->size(); n++) {
        particle_t& particle = allParticles->at(n);
        idPosPairs.push_back(
            std::make_tuple(std::addressof(particle), particle.position[ndim]));
    }
    std::sort(idPosPairs.begin(), idPosPairs.end(), _CompParticlesPos());
    int groupSize = allParticles->size() / ngroup;
    groupParticles->resize(ngroup);
    int np = 0;
    for(int ng=0; ng<ngroup; ng++) {
        for(int ngp=0; ngp<groupSize && np<allParticles->size(); ngp++) {
            groupParticles->at(ng).push_back(*std::get<0>(idPosPairs.at(np)));
            np++;
        }
    }
    for(; np<allParticles->size(); np++) {
        groupParticles->at(ngroup-1).push_back(*std::get<0>(idPosPairs.at(np)));
    }
}

void dispatch(std::vector<particle_t>* allParticles, 
                std::vector<std::vector<particle_t>>* localParticleGroups,
                topology_t* topology
            ) {
    int& gridX = topology->gridX;
    int& gridY = topology->gridY;
    int& gridZ = topology->gridZ;
    int& subgridX = topology->subgridX;
    int& subgridY = topology->subgridY;
    int& subgridZ = topology->subgridZ;
    int& rankX = topology->rankX;
    int& rankY = topology->rankY;
    int& rankZ = topology->rankZ;

    int totalSX = gridX*subgridX;
    std::vector<std::vector<particle_t>> slices;
    std::vector<particle_t> slice;
    std::vector<std::vector<particle_t>> fibers;
    std::vector<particle_t> fiber;
    std::vector<std::vector<particle_t>> blocks;
    std::vector<particle_t> block;
    if(rankY==0 && rankZ==0) {
        if(rankX==0) {
            slices.resize(gridX);
            for(int i=0; i<gridX; i++)
                slices.at(i).clear();
            divideGroupN(allParticles, &slices, 0, gridX);
            for(int i=0; i<slices.size(); i++) {
                if(i==0) {
                    slice.clear();
                    for(int j=0; j<slices.at(i).size(); j++)
                        slice.push_back(slices.at(i).at(j));
                }
                else {
                    int sizeSlice = slices.at(i).size();
                    MPI_Send(
                        &sizeSlice,
                        1,
                        MPI_INT,
                        i,
                        i,
                        topology->xEdgeComm
                    );
                    MPI_Send(
                        const_cast<particle_t*>(slices.at(i).data()),
                        sizeSlice*sizeof(particle_t),
                        MPI_BYTE,
                        i,
                        i,
                        topology->xEdgeComm
                    );
                }
            }
        }
        else {
            int sizeSlice = 0;
            MPI_Recv(&sizeSlice, 1, MPI_INT, 0, rankX, topology->xEdgeComm, nullptr);
            slice.resize(sizeSlice);
            MPI_Recv(
                const_cast<particle_t*>(slice.data()), 
                sizeSlice*sizeof(particle_t),
                MPI_BYTE,
                0, 
                rankX,
                topology->xEdgeComm,
                nullptr
            );
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(rankZ==0) {
        if(rankY==0) {
            fibers.resize(gridY);
            for(int i=0; i<gridY; i++)
                fibers.at(i).clear();
            divideGroupN(&slice, &fibers, 1, gridY);
            for(int i=0; i<fibers.size(); i++) {
                if(i==0) {
                    fiber.clear();
                    for(int j=0; j<fibers.at(i).size(); j++)
                        fiber.push_back(fibers.at(i).at(j));
                }
                else {
                    int sizeFiber = fibers.at(i).size();
                    MPI_Send(
                        &sizeFiber,
                        1,
                        MPI_INT,
                        i,
                        i,
                        topology->yEdgeComm
                    );
                    MPI_Send(
                        const_cast<particle_t*>(fibers.at(i).data()),
                        sizeFiber*sizeof(particle_t),
                        MPI_BYTE,
                        i,
                        i,
                        topology->yEdgeComm
                    );
                }
            }
        }
        else {
            int sizeFiber = 0;
            MPI_Recv(&sizeFiber, 1, MPI_INT, 0, rankY, topology->yEdgeComm, nullptr);
            fiber.resize(sizeFiber);
            MPI_Recv(
                const_cast<particle_t*>(fiber.data()),
                sizeFiber*sizeof(particle_t),
                MPI_BYTE,
                0,
                rankY,
                topology->yEdgeComm,
                nullptr
            );
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(rankZ==0) {
        blocks.resize(gridZ);
        for(int i=0; i<gridZ; i++) 
            blocks.at(i).clear();
        divideGroupN(&fiber, &blocks, 2, gridZ);
        for(int i=0; i<blocks.size(); i++) {
            if(i==0) {
                block.clear();
                for(int j=0; j<blocks.at(i).size(); j++)
                    block.push_back(blocks.at(i).at(j));
            }
            else {
                int sizeBlock = blocks.at(i).size();
                MPI_Send(
                    &sizeBlock,
                    1,
                    MPI_INT,
                    i,
                    i,
                    topology->zEdgeComm
                );
                MPI_Send(
                    const_cast<particle_t*>(blocks.at(i).data()),
                    sizeBlock*sizeof(particle_t),
                    MPI_BYTE,
                    i,
                    i,
                    topology->zEdgeComm
                );
            }
        }
    }
    else {
        int sizeBlock = 0;
        MPI_Recv(&sizeBlock, 1, MPI_INT, 0, rankZ, topology->zEdgeComm, nullptr);
        block.resize(sizeBlock);
        MPI_Recv(
            const_cast<particle_t*>(block.data()),
            sizeBlock*sizeof(particle_t),
            MPI_BYTE,
            0,
            rankZ,
            topology->zEdgeComm,
            nullptr
        );
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // block.clear();
    for(int i=0; i<blocks.size(); i++)
        blocks.at(i).clear();
    blocks.clear();

    fiber.clear();
    for(int i=0; i<fibers.size(); i++)
        fibers.at(i).clear();
    fibers.clear();

    slice.clear();
    for(int i=0; i<slices.size(); i++)
        slices.at(i).clear();
    slices.clear();
    localParticleGroups->resize(0);


    divideGroupN(&block, &slices, 0, subgridX);
    for(int sgx=0; sgx<subgridX; sgx++) {
        divideGroupN(std::addressof(slices.at(sgx)), &fibers, 1, subgridY);
        for(int sgy=0; sgy<subgridY; sgy++) {
            divideGroupN(std::addressof(fibers.at(sgy)), &blocks, 2, subgridZ);
            for(int sgz=0; sgz<subgridZ; sgz++) {
                localParticleGroups->push_back(blocks.at(sgz));
            }
        }
    }
    
    block.clear();
    for(int i=0; i<blocks.size(); i++)
        blocks.at(i).clear();
    blocks.clear();

    fiber.clear();
    for(int i=0; i<fibers.size(); i++)
        fibers.at(i).clear();
    fibers.clear();

    slice.clear();
    for(int i=0; i<slices.size(); i++)
        slices.at(i).clear();
    slices.clear();

    std::vector<particle_t> blockCentroidsTx;
    std::vector<particle_t> blockCentroidsRx;
    blockCentroidsTx.resize(gridX*gridY*gridZ*subgridX*subgridY*subgridZ);
    blockCentroidsRx.resize(gridX*gridY*gridZ*subgridX*subgridY*subgridZ);
    for(int nn=0; nn<subgridX*subgridY*subgridZ; nn++) {
        particle_t* p = std::addressof(blockCentroidsTx.at(nn));
        getCentroid(std::addressof(localParticleGroups->at(nn)), &p, 0);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(
        const_cast<particle_t*>(blockCentroidsTx.data()),
        subgridX*subgridY*subgridZ*sizeof(particle_t),
        MPI_BYTE,
        const_cast<particle_t*>(blockCentroidsRx.data()),
        subgridX*subgridY*subgridZ*sizeof(particle_t),
        MPI_BYTE,
        MPI_COMM_WORLD
    );
    for(int n=0; n<blockCentroidsRx.size(); n++) {
        for(int nn=0; nn<localParticleGroups->size(); nn++)
            localParticleGroups->at(nn).push_back(blockCentroidsRx.at(n));
    }
    blockCentroidsRx.clear();
    blockCentroidsTx.clear();
}


void gather(std::vector<std::vector<particle_t>>* localParticleGroups, 
            std::vector<particle_t>* allParticles, 
            topology_t* topology) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(rank==0) {
        allParticles->clear();
        for(int i=0; i<size; i++) {
            if(i==0) {
                for(int n=0; n<localParticleGroups->size(); n++)
                    for(int nn=0; nn<localParticleGroups->at(n).size(); nn++)
                        allParticles->push_back(localParticleGroups->at(n).at(nn));
            }
            else {
                std::vector<particle_t> particlesRx;
                int ngroup;
                MPI_Recv(&ngroup, 1, MPI_INT, i, 0, MPI_COMM_WORLD, nullptr);
                for(int n=0; n<ngroup; n++) {
                    int groupSize;
                    MPI_Recv(&groupSize, 1, MPI_INT, i, 2*n+1, MPI_COMM_WORLD, nullptr);
                    particlesRx.resize(groupSize);
                    MPI_Recv(
                        const_cast<particle_t*>(particlesRx.data()),
                        groupSize * sizeof(particle_t),
                        MPI_BYTE,
                        i, 
                        2*n+2,
                        MPI_COMM_WORLD,
                        nullptr
                    );
                    for(int j=0; j<particlesRx.size(); j++){
                        if(particlesRx.at(j).id&0xff00000000000000ul)
                            continue;
                        allParticles->push_back(particlesRx.at(j));
                    }
                    particlesRx.clear();
                }
            }
        }
    }
    else {
        int ngroups = localParticleGroups->size();
        MPI_Send(
            &ngroups,
            1,
            MPI_INT,
            0,
            0,
            MPI_COMM_WORLD
        );
        for(int n=0; n<ngroups; n++) {
            int groupSize = localParticleGroups->at(n).size();
            MPI_Send(
                &groupSize,
                1,
                MPI_INT,
                0,
                2*n+1,
                MPI_COMM_WORLD
            );
            MPI_Send(
                const_cast<particle_t*>(localParticleGroups->at(n).data()),
                groupSize*sizeof(particle_t),
                MPI_BYTE,
                0,
                2*n+2,
                MPI_COMM_WORLD
            );

        }
    }
    return;
}

void getCentroid(std::vector<particle_t>* particles, particle_t** centroid, int weightFeat) {
    uint64_t id = 0;
    if(*centroid == nullptr)
        *centroid = static_cast<particle_t*>(malloc(sizeof(particle_t)));
    assert(weightFeat < PARTICLE_N_FEAT);
    double pos[PARTICLE_N_DIM];
    double vel[PARTICLE_N_DIM];
    double feat[PARTICLE_N_FEAT];
    for(int n=0; n<PARTICLE_N_DIM; n++) {
        pos[n] = 0.0;
        vel[n] = 0.0;
    }
    for(int n=0; n<PARTICLE_N_FEAT; n++)
        feat[n] = 0.0;
    for(int n=0; n<particles->size(); n++) {
        particle_t& particle = particles->at(n);
        for(int i=0; i<PARTICLE_N_DIM; i++) {
            pos[i] += particle.position[i] * particle.features[weightFeat];
            vel[i] += particle.velocity[i] * particle.features[weightFeat];
        }
        for(int i=0; i<PARTICLE_N_FEAT; i++)
            feat[i] += particle.features[i];
        id += particle.id;
    }
    for(int n=0; n<PARTICLE_N_DIM; n++) {
        pos[n] /= feat[weightFeat];
        vel[n] /= feat[weightFeat];
        (*centroid)->position[n] = pos[n];
        (*centroid)->velocity[n] = vel[n];
        (*centroid)->acceleration[n] = 0.0;
    }
    for(int n=0; n<PARTICLE_N_FEAT; n++) {
        (*centroid)->features[n] = feat[n];
    }
    (*centroid)->enabled = true;
    (*centroid)->id = (id>>8) | 0xff00000000000000ul;
    
}
