#include "particle.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <cassert>

void updateVelocity(particle_t* particle, double timeStep) {
    double* accelerationDim = particle->acceleration;
    double* velocityDim = particle->velocity;
    for(int n=0; n<particle->ndim; n++) 
        *(velocityDim++) += *(accelerationDim++) * timeStep;
}

void updatePosition(particle_t* particle, double timeStep) {
    double* positionDim = particle->position;
    double* velocityDim = particle->velocity;
    for(int n=0; n<particle->ndim; n++) 
        *(positionDim++) += *(velocityDim++) * timeStep;
}

void updateAcceleration(particle_t* one, particle_t* another) {
    (*(one->updateAcceleration))(one, another);
}

void load_particles(const std::string filename, particle_t** particles, int* nParticles) {
    std::ifstream file(filename, std::ios::in);
    if(!file.is_open()) {
        std::fprintf(stderr, "[%s:%d]: can not open the file %s\n", __FILE__, __LINE__, filename);
        exit(EXIT_FAILURE);
    }
    int space_dimension = -1;
    std::string line;
    std::getline(file, line);
    std::sscanf(line.c_str(), "%d %d", nParticles, &space_dimension);
    printf("%d\n", *nParticles);
    particle_t* particles_ = static_cast<particle_t*>(malloc(*nParticles * sizeof(particle_t)));
    double* pos = static_cast<double*>(malloc(*nParticles * space_dimension * sizeof(double)));
    assert(pos != nullptr);
    double* vel = static_cast<double*>(malloc(*nParticles * space_dimension * sizeof(double)));
    assert(vel != nullptr);
    double* acc = static_cast<double*>(malloc(*nParticles * space_dimension * sizeof(double)));
    assert(acc != nullptr);
    double* features = nullptr;
    for(int n=0; n<*nParticles; n++) {
        (particles_+n)->position = pos+n*space_dimension;
        (particles_+n)->velocity = vel+n*space_dimension;
        (particles_+n)->acceleration = acc+n*space_dimension;
        (particles_+n)->ndim = space_dimension;
        (particles_+n)->updateAcceleration = nullptr;
        (particles_+n)->features = nullptr;
        (particles_+n)->nfeat = -1;
    }
    particle_t* pp = particles_;
    std::stringstream lineStream;
    double foo;
    int n = 0;
    int nfeat = -1;
    while(std::getline(file, line) && n<*nParticles) {
        printf("reading: %d\n", n);
        lineStream.str(line);
        lineStream >> foo;
        pp->id = static_cast<uint64_t>(foo);
        for(int nn=0; nn<pp->ndim; nn++)
            lineStream >> *(pp->position+nn);
        for(int nn=0; nn<pp->ndim; nn++)
            lineStream >> *(pp->velocity+nn);
        for(int nn=0; nn<pp->ndim; nn++)
            lineStream >> *(pp->acceleration+nn);
        lineStream >> foo;
        pp->nfeat = static_cast<int>(foo);
        // lineStream >> pp->nfeat;
        if (nfeat == -1) {
            nfeat = pp->nfeat;
            features = static_cast<double*>(malloc(nfeat*(*nParticles)*sizeof(double)));
        }
        else {
            assert(nfeat==pp->nfeat);
        }
        pp->features = features + n*nfeat;
        for(int nn=0; nn<pp->nfeat; nn++)
            lineStream >> *(pp->features+nn);
        n++;
    }
    file.close();
    *particles = particles_;
    if(n != *nParticles) {
        std::fprintf(stderr, "[%s:%d]: unmatched particle number %d vs %d\n", __FILE__, __LINE__, n, *nParticles);
        exit(EXIT_FAILURE);
    }
}

void store_particles(const std::string filename, const particle_t* particles, const int nParticles) {
    std::ofstream file(filename, std::ios::out|std::ios::trunc);
    if(!file.is_open()) {
        std::fprintf(stderr, "[%s:%d]: can not open the file %s\n", __FILE__, __LINE__, filename);
        exit(EXIT_FAILURE);
    }
    int space_dimension = particles->ndim;
    file << nParticles << " " << space_dimension << std::endl;
    file.flush();
    particle_t* pp = const_cast<particle_t*>(particles);
    for(int n=0; n<nParticles; n++) {
        printf("storing: %d\n", n);
        file << pp->id << " ";
        for(int nn=0; nn<pp->ndim; pp++)
            file << *((pp->position)+nn) << " ";
        for(int nn=0; nn<pp->ndim; pp++)
            file << *((pp->velocity)+nn) << " ";
        for(int nn=0; nn<pp->ndim; pp++)
            file << *((pp->acceleration)+nn) << " ";
        file << pp->nfeat << " ";
        for(int nn=0; nn<pp->nfeat; nn++) {
            file << *(pp->features+nn);
            if(nn+1 != pp->nfeat)
                file << " ";
            else
                file << std::endl;
        }
        pp++;
    }
    file.close();
}

