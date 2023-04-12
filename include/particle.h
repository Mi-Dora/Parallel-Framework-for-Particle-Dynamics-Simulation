#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include <cstdint>
#include <cstdlib>

struct particle_t {
    double* position;
    double* velocity;
    double* acceleration;
    double* features;
    void (*updateAcceleration)(particle_t*, particle_t*);
    int ndim;
};

void updateVelocity(particle_t* particle, double timeStep);
void updatePosition(particle_t* particle, double timeStep);
void updateAcceleration(particle_t* one, particle_t* another);

#endif
