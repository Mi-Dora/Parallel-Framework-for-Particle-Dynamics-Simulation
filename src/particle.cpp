#include "particle.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

void updateVelocity(particle_t* particle, double timeStep) {
    double* accelerationDim = particle->acceleration;
    double* velocityDim = particle->velocity;
    for(size_t n=0; n<particle->ndim; n++) 
        *(velocityDim++) += *(accelerationDim++) * timeStep;
}

void updatePosition(particle_t* particle, double timeStep) {
    double* positionDim = particle->position;
    double* velocityDim = particle->velocity;
    for(size_t n=0; n<particle->ndim; n++) 
        *(positionDim++) += *(velocityDim++) * timeStep;
}

void updateAcceleration(particle_t* one, particle_t* another) {
    (*(one->updateAcceleration))(one, another);
}

void load_particles(const std::string filename, particle_t** particles, size_t* nParticles) {
    std::ifstream file(filename, std::ios::in);
    if(!file.is_open()) {
        std::fprintf(stderr, "[%s:%d]: can not open the file %s\n", __FILE__, __LINE__, filename);
        exit(EXIT_FAILURE);
    }
    int space_dimension = -1;
    std::string line;
    std::getline(file, line);
    std::sscanf(line.c_str(), "%d %d", nParticles, &space_dimension);
    particle_t* particles_ = static_cast<particle_t*>(malloc(*nParticles * sizeof(particle_t)));
    double* pos = static_cast<double*>(malloc(*nParticles * space_dimension * sizeof(double)));
    double* vel = static_cast<double*>(malloc(*nParticles * space_dimension * sizeof(double)));
    double* acc = static_cast<double*>(malloc(*nParticles * space_dimension * sizeof(double)));
    for(size_t n=0; n<*nParticles; n++) {
        (particles_+n)->position = pos+n*space_dimension*sizeof(double);
        (particles_+n)->velocity = vel+n*space_dimension*sizeof(double);
        (particles_+n)->acceleration = acc+n*space_dimension*sizeof(double);
        (particles_+n)->ndim = space_dimension;
        (particles_+n)->updateAcceleration = nullptr;
        (particles_+n)->features = nullptr;
        (particles_+n)->nfeat = -1;
    }
    particle_t* pp = particles_;
    std::stringstream lineStream;
    size_t n = 0;
    while(std::getline(file, line) && n<*nParticles) {
        n++;
        lineStream.str(line);
        for(size_t nn=0; nn<pp->ndim; nn++)
            lineStream >> *(pp->position+nn);
        for(size_t nn=0; nn<pp->ndim; nn++)
            lineStream >> *(pp->velocity+nn);
        for(size_t nn=0; nn<pp->ndim; nn++)
            lineStream >> *(pp->acceleration+nn);
        // int feature_dimension = -1;
        lineStream >> pp->nfeat;
        pp->features = static_cast<double*>(malloc(pp->nfeat*sizeof(double)));
        for(size_t nn=0; nn<pp->nfeat; nn++)
            lineStream >> *(pp->features+nn);
    }
    if(n != *nParticles) {
        std::fprintf(stderr, "[%s:%d]: unmatched particle number %d vs %d\n", __FILE__, __LINE__, n, *nParticles);
        exit(EXIT_FAILURE);
    }
}

