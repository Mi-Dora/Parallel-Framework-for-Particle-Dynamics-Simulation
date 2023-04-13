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

struct chunk_particles_t {
    particle_t* particles;
    int nParticle = 0;  
};

void updateVelocity(particle_t* particle, double timeStep);
void updatePosition(particle_t* particle, double timeStep);
void updateAcceleration(particle_t* one, particle_t* another);
void clearAccelerations(particle_t* particles, int nParticle);
void updates(particle_t* particle, int nParticle, double timeStep);
void registerUpdateAccelerationFn(particle_t* one, void (*fn)(particle_t*, particle_t*));

void load_particles(const std::string filename, particle_t** particles, int* nParticles);
void load_particles(const std::string filename, chunk_particles_t** particleChunk);
void store_particles(const std::string filename, particle_t* particles, int nParticles);
void store_particles(const std::string filename, chunk_particles_t* particleChunk);

void free_particles(particle_t* particles);
void free_particles(chunk_particles_t* particles);

void updates(particle_t* particles, int nParticle, double timeStep);
void updates(chunk_particles_t* particlesChunk, double timeStep);

#endif
