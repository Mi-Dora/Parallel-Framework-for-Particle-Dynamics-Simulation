#ifndef __NEO_PARTICLE_H__
#define __NEO_PARTICLE_H__

#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

#define PARTICLE_N_DIM 3
#define PARTICLE_N_FEAT 5

struct particle_t {
    uint64_t id;
    double position[PARTICLE_N_DIM];
    double velocity[PARTICLE_N_DIM];
    double acceleration[PARTICLE_N_DIM];
    double features[PARTICLE_N_FEAT];
    bool enabled;
};

void updateVelocity(particle_t* particle, double timeStep);
void updatePosition(particle_t* particle, double timeStep);
void clearAcceleration(particle_t* particle);

void loadParticles(const std::string filename, std::vector<particle_t>* particles);
void storeParticles(const std::string filename, std::vector<particle_t>* particles);

#endif
