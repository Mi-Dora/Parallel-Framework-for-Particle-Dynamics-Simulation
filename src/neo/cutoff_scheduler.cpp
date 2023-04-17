#include "neo/cutoff_scheduler.h"
#include "neo/particle.h"
#include <vector>
#include <algorithm>
#include <tuple>

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
                std::vector<particle_t>* localParticles,
                std::vector<particle_t>* neighborParticles,
                topology_t* topology
                ) {
    int& gridX = topology->gridX;
    int& gridY = topology->gridY;
    int& gridZ = topology->gridZ;
    int& subgridX = topology->subgridX;
    int& subgridY = topology->subgridY;
    int& subgridZ = topology->subgridZ;
    
}
