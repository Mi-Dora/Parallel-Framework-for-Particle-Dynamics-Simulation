#include "neo/particle.h"
#include "neo/dense_scheduler.h"
#include "neo/user_def/gravity_particle.h"
#include <string>
#include "mpi.h"
#include <vector>

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
    std::vector<particle_t> globalParticles;
    
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
    if(rankX==0 && rankY==0) {
        loadParticles(inputFilename, &globalParticles);
        padding(&globalParticles, minimalPaddingSize(globalParticles.size(), gridX, gridY, 1));
        storeParticles(paddingFilename, &globalParticles);
    }

    std::vector<particle_t> localParticlesX, localParticlesY;
    scatter(&localParticlesX, &localParticlesY, &globalParticles, &topology);
    for(int iter=0; iter<100; iter++) {
        update(&localParticlesX, &localParticlesY, &topology, 1e-1, gravityUpdateAcceleration);
    }

    storeParticles(scatterFilenameX, &localParticlesX);
    storeParticles(scatterFilenameY, &localParticlesY);
    gather(&localParticlesX, &localParticlesY, &globalParticles, &topology);
    if(rankX==0 && rankY==0){
        storeParticles(gatherFilename, &globalParticles);
        shrink(&globalParticles);
        storeParticles(shrinkFilename, &globalParticles);
    }

    MPI_Finalize();
    return 0;
}
