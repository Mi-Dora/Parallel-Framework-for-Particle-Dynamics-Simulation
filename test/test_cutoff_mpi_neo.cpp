#include "neo/particle.h"
#include "neo/cutoff_scheduler.h"
#include "neo/user_def/gravity_particle.h"
#include <string>
#include "mpi.h"
#include <vector>
#include <cassert>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int gridX, gridY, gridZ;
    assert(argc==5);
    gridX = atoi(argv[1]);
    gridY = atoi(argv[2]);
    gridZ = atoi(argv[3]);
    std::string inputFilename(argv[4]);

    topology_t topology;
    setupCommunicators(&topology, gridX, gridY, gridZ, 1, 1, 1);
    const int& rankX = topology.rankX;
    const int& rankY = topology.rankY;
    const int& rankZ = topology.rankZ;
    std::vector<particle_t> globalParticles;
    
    if(rankX==0 && rankY==0) {
        loadParticles(inputFilename, &globalParticles);
        // padding(&globalParticles, minimalPaddingSize(globalParticles.size(), gridX, gridY, 1));
        // storeParticles(paddingFilename, &globalParticles);
    }
    std::vector<std::vector<particle_t>> localParticleGroups;
    dispatch(&globalParticles, &localParticleGroups, &topology);

    double time = 0, totalTime = 0;
    totalTime = -MPI_Wtime();
    constexpr int totalIter = 100;
    for(int iter=0; iter<totalIter; iter++) {
        for(int n=0; n<localParticleGroups.size(); n++) {
            particle_t* particles = const_cast<particle_t*>(localParticleGroups.at(n).data());
            for(int nx=0; nx<localParticleGroups.at(n).size(); nx++)
                clearAcceleration(
                    particles + nx
                );
            for(int nx=0; nx<localParticleGroups.at(n).size(); nx++) {
                particle_t* pa = std::addressof(localParticleGroups.at(n).at(nx));
                for(int ny=0; ny<localParticleGroups.at(n).size(); ny++) {
                    particle_t* pb = std::addressof(localParticleGroups.at(n).at(ny));
                    if(pa->id == pb->id)
                        continue;
                    gravityUpdateAcceleration(pa, pb);
                }
            }
            for(int nx=0; nx<localParticleGroups.at(n).size(); nx++){
                updateVelocity(particles+nx, 1e-1);
                updatePosition(particles+nx, 1e-1);
            }
        }
        // if(rankX==0 && rankY==0) {
        //     time = -MPI_Wtime();
        //     printf("iter %d", iter);
        // }
        // update(&localParticlesX, &localParticlesY, &topology, 1e-1, gravityUpdateAcceleration);
        // if(rankX==0 && rankY==0) {
        //     time += MPI_Wtime();
        //     printf("iter time %f\n", time);
        // }
    }
    totalTime += MPI_Wtime();
    if(rankX==0 && rankY==0 && rankZ==0)
        printf("avg iter time %f\n", totalTime / totalIter);
    gather(&localParticleGroups, &globalParticles, &topology);

    // storeParticles(scatterFilenameX, &localParticlesX);
    // storeParticles(scatterFilenameY, &localParticlesY);
    if(rankX==0 && rankY==0){
        // storeParticles(gatherFilename, &globalParticles);
        // shrink(&globalParticles);
        // storeParticles(shrinkFilename, &globalParticles);
    }

    MPI_Finalize();
    return 0;
}
