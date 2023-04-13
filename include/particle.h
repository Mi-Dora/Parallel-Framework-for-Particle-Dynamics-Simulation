#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include <cstdint>
#include <cstdlib>
#include <string>

struct particle_t {
    uint64_t id;
    double* position;
    double* velocity;
    double* acceleration;
    double* features;
    void (*updateAcceleration)(particle_t*, particle_t*);
    int ndim;
    int nfeat;
};

void updateVelocity(particle_t* particle, double timeStep);
void updatePosition(particle_t* particle, double timeStep);
void updateAcceleration(particle_t* one, particle_t* another);

void load_particles(const std::string filename, particle_t** particles, int* nParticles);
void store_particles(const std::string filename, const particle_t* particles, const int nParticles);

#endif
