#include "particle.h"
#include "dense_scheduler.h"
#include "user_def/gravity_particle.h"
#include <string>
#include "mpi.h"

int main(int argc, char** argv) {
    const int gridX = 3, gridY=3;
    MPI_Init(&argc, &argv);

    topology_t topology;
    setupCommunicators(&topology, gridX, gridY);
    const int& rankX = topology.rankX;
    const int& rankY = topology.rankY;
    const std::string inputFilename = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/data/data100.txt";
    const std::string paddingFilename = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/data/padding100.txt";
    const std::string scatterFilenameX = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/data/scatter100.txt.x" + std::to_string(rankX) + std::to_string(rankY);
    const std::string scatterFilenameY = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/data/scatter100.txt.y" + std::to_string(rankX) + std::to_string(rankY);
    const std::string gatherFilename = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/data/gather100.txt";
    const std::string shrinkFilename = "/home/cman8/Parallel-Framework-for-Particle-Dynamics-Simulation/data/shrink100.txt";
    chunk_particles_t* globalParticles = nullptr;
    
    int sizeReduceXComm=-1, sizeReduceYComm=-1, sizeScatterXComm=-1, sizeScatterYComm=-1;
    if(topology.reduceXComm != MPI_COMM_NULL)
        MPI_Comm_size(topology.reduceXComm, &sizeReduceXComm);
    if(topology.reduceYComm != MPI_COMM_NULL)
        MPI_Comm_size(topology.reduceYComm, &sizeReduceYComm);
    if(topology.scatterXComm != MPI_COMM_NULL)
        MPI_Comm_size(topology.scatterXComm, &sizeScatterXComm);
    if(topology.scatterYComm != MPI_COMM_NULL)
        MPI_Comm_size(topology.scatterYComm, &sizeScatterYComm);
    std::printf(
        "(%d, %d) => (%d, %d) (%lx:%d %lx:%d %lx:%d %lx:%d)\n", 
        topology.gridX, topology.gridY, topology.rankX, topology.rankY, 
        topology.reduceXComm, sizeReduceXComm, topology.reduceYComm, sizeReduceYComm,
        topology.scatterXComm, sizeScatterXComm, topology.scatterYComm, sizeScatterYComm);
    if(globalParticles!=nullptr)
        printf("FSFS before_gather==%lx\n", globalParticles->particles);
    else
        printf("FSFS before_gather=null\n");
    if(rankX==0 && rankY==0) {
        load_particles(inputFilename, &globalParticles);
        padding(globalParticles, gridX*gridY);
        store_particles(paddingFilename, globalParticles);
    }
    if(globalParticles!=nullptr)
        printf("FSKS before_gather==%lx\n", globalParticles->particles);
    else
        printf("FSKS before_gather=null\n");
    chunk_particles_t *localParticlesX=nullptr, *localParticlesY=nullptr;
    scatter(&localParticlesX, &localParticlesY, globalParticles, &topology);
    store_particles(scatterFilenameX, localParticlesX);
    store_particles(scatterFilenameY, localParticlesY);
    if(globalParticles!=nullptr)
        printf("before_gather==%lx", globalParticles->particles);
    else
        printf("before_gather=null");
    gather(localParticlesX, localParticlesY, globalParticles, &topology);
    if(globalParticles!=nullptr)
        printf(" after_gather==%lx\n", globalParticles->particles);
    else
        printf(" after_gather==nullptr\n");
    if(rankX==0 && rankY==0){
        store_particles(gatherFilename, globalParticles);
        shrink(globalParticles);
        store_particles(shrinkFilename, globalParticles);
    }

    MPI_Finalize();
    return 0;
}
