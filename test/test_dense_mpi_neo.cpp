#include "neo/particle.h"
#include "neo/dense_scheduler.h"
#include "neo/user_def/gravity_particle.h"
#include <string>
#include "mpi.h"
#include <vector>
#include <cassert>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int gridX, gridY;
    assert(argc==4);
    gridX = atoi(argv[1]);
    gridY = atoi(argv[2]);
    std::string inputFilename(argv[3]);

    topology_t topology;
    setupCommunicators(&topology, gridX, gridY);
    const int& rankX = topology.rankX;
    const int& rankY = topology.rankY;
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
        // storeParticles(paddingFilename, &globalParticles);
    }

    std::vector<particle_t> localParticlesX, localParticlesY;
    scatter(&localParticlesX, &localParticlesY, &globalParticles, &topology);
    double time = 0, totalTime = 0;
    totalTime = -MPI_Wtime();
    constexpr int totalIter = 100;
    for(int iter=0; iter<totalIter; iter++) {
        // if(rankX==0 && rankY==0) {
        //     time = -MPI_Wtime();
        //     printf("iter %d", iter);
        // }
        update(&localParticlesX, &localParticlesY, &topology, 1e-1, gravityUpdateAcceleration);
        // if(rankX==0 && rankY==0) {
        //     time += MPI_Wtime();
        //     printf("iter time %f\n", time);
        // }
    }
    totalTime += MPI_Wtime();
    printf("avg iter time %f\n", totalTime / totalIter);

    // storeParticles(scatterFilenameX, &localParticlesX);
    // storeParticles(scatterFilenameY, &localParticlesY);
    gather(&localParticlesX, &localParticlesY, &globalParticles, &topology);
    if(rankX==0 && rankY==0){
        // storeParticles(gatherFilename, &globalParticles);
        shrink(&globalParticles);
        // storeParticles(shrinkFilename, &globalParticles);
    }

    MPI_Finalize();
    return 0;
}
