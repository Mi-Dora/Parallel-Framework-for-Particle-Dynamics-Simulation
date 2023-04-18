#include "neo/cutoff_scheduler.h"
#include "neo/particle.h"
#include <vector>
#include <algorithm>
#include <tuple>
#include <mpi.h>

class _CompParticlesPos {
public:
    bool operator() (std::tuple<particle_t*, double> a, std::tuple<particle_t*, double> b) {
        return (a<b);
    }
};

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
                    MPI_INIT,
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

    // TODO: no subgrid, need to complete
    localParticleGroups->resize(0);
    localParticleGroups->at(0).push_back(block);

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
}


void gather(std::vector<std::vector<particle_t>>* localParticleGroups, 
            std::vector<particle_t>* allParticles, 
            topology_t* topology) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank==0) {

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
}
