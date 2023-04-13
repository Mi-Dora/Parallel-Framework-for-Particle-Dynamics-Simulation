#include "particle.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <cassert>
#include "string.h"
#include "stdio.h"

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
    if (one->updateAcceleration == nullptr) 
        std::cerr << "no updateAcceleration registerd on particle: " << one->id << std::endl;
    else
        (*(one->updateAcceleration))(one, another);
}

void registerUpdateAccelerationFn(particle_t* one, void (*fn)(particle_t*, particle_t*)) {
    one->updateAcceleration = fn;
}

void registerUpdateAccelerationFns(particle_t* particles, int nParticles, void (*fn)(particle_t*, particle_t*)) {
    particle_t* pt = particles;
    for (int n=0; n<nParticles; n++)
        registerUpdateAccelerationFn(pt++, fn);
}

void registerUpdateAccelerationFns(chunk_particles_t* particleChunk, void (*fn)(particle_t*, particle_t*)) {
    registerUpdateAccelerationFns(particleChunk->particles, particleChunk->nParticle, fn);
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
    particle_t* particles_ = static_cast<particle_t*>(malloc(*nParticles * sizeof(particle_t)));
    double* pos = static_cast<double*>(malloc((*nParticles) * space_dimension * sizeof(double)));
    assert(pos != nullptr);
    double* vel = static_cast<double*>(malloc((*nParticles) * space_dimension * sizeof(double)));
    assert(vel != nullptr);
    double* acc = static_cast<double*>(malloc((*nParticles) * space_dimension * sizeof(double)));
    assert(acc != nullptr);
    double* features = nullptr;
    particle_t* pp = particles_;
    for(int n=0; n<(*nParticles); n++) {
        pp->position = pos+n*space_dimension;
        pp->velocity = vel+n*space_dimension;
        pp->acceleration = acc+n*space_dimension;
        pp->ndim = space_dimension;
        pp->updateAcceleration = nullptr;
        pp->features = nullptr;
        pp->nfeat = -1;
        pp++;
    }

    pp = particles_;
    std::stringstream lineStream;
    double foo;
    int n = 0;
    int nfeat = -1;
    while(std::getline(file, line) && n<(*nParticles)) {
        char* bar;
        bar = strtok(const_cast<char*>(line.c_str()), " ");
        foo = atof(bar);
        pp->id = static_cast<uint64_t>(foo);
        for(int nn=0; nn<(pp->ndim); nn++) {
            bar = strtok(NULL, " ");
            *(pp->position + nn) = atof(bar);
        }
        for(int nn=0; nn<(pp->ndim); nn++) {
            bar = strtok(NULL, " ");
            *(pp->velocity + nn) = atof(bar);
        }
        for(int nn=0; nn<(pp->ndim); nn++) {
            bar = strtok(NULL, " ");
            *(pp->acceleration + nn) = atof(bar);
        }
        bar = strtok(NULL, " ");
        foo = atof(bar);
        pp->nfeat = static_cast<int>(foo);
        if (nfeat == -1) {
            nfeat = pp->nfeat;
            features = static_cast<double*>(malloc(nfeat*(*nParticles)*sizeof(double)));
        }
        else {
            assert(nfeat==(pp->nfeat));
        }
        pp->features = features + n*nfeat;
        for(int nn=0; nn<(pp->nfeat); nn++) {
            bar = strtok(NULL, " ");
            *(pp->features + nn) = atof(bar);
        }
        n++;
        pp++;
    }
    file.close();
    *particles = particles_;
    if(n != *nParticles) {
        std::fprintf(stderr, "[%s:%d]: unmatched particle number %d vs %d\n", __FILE__, __LINE__, n, *nParticles);
        exit(EXIT_FAILURE);
    }
}

void store_particles(const std::string filename, particle_t* particles, int nParticles) {
    std::ofstream file(filename, std::ios::out|std::ios::trunc);
    if(!file.is_open()) {
        std::fprintf(stderr, "[%s:%d]: can not open the file %s\n", __FILE__, __LINE__, filename);
        exit(EXIT_FAILURE);
    }
    int space_dimension = particles->ndim;
    file << nParticles << " " << space_dimension << std::endl;
    file.flush();
    particle_t* pp = particles;
    for(int n=0; n<nParticles; n++) {
        pp = particles + n;
        file << std::scientific << pp->id << " ";
        for(int nn=0; nn<(pp->ndim); nn++)
            file << std::scientific << *((pp->position)+nn) << " ";
        for(int nn=0; nn<(pp->ndim); nn++)
            file << std::scientific << *((pp->velocity)+nn) << " ";
        for(int nn=0; nn<(pp->ndim); nn++)
            file << std::scientific << *((pp->acceleration)+nn) << " ";
        file << pp->nfeat << " ";
        for(int nn=0; nn<pp->nfeat; nn++) {
            file << std::scientific << *(pp->features+nn);
            if(nn+1 != pp->nfeat)
                file << " ";
            else
                file << std::endl;
        }
    }
    file.close();
}

void load_particles(const std::string filename, chunk_particles_t** particleChunk) {
    if(*particleChunk == nullptr)
        *particleChunk = static_cast<chunk_particles_t*>(malloc(sizeof(particleChunk)));
    else
        assert((*particleChunk)->nParticle == 0);
    particle_t* particles;
    int nParticles;
    load_particles(filename, &particles, &nParticles);
    (*particleChunk)->particles = particles;
    (*particleChunk)->nParticle = nParticles;
}

void store_particles(const std::string filename, chunk_particles_t* particleChunk) {
    store_particles(filename, particleChunk->particles, particleChunk->nParticle);
}

void free_particles(particle_t* particles) {
    free(particles->position);
    free(particles->velocity);
    free(particles->acceleration);
    free(particles->features);
    free(particles);
}

void free_particles(chunk_particles_t* particles) {
    free(particles->particles);
    particles->particles = nullptr;
    particles->nParticle = 0;
}

void clearAccelerations(particle_t* particles, int nParticle) {
    particle_t* pt = particles;
    for(int n=0; n<nParticle; n++) {
        double* acceleration = pt->acceleration;
        for(int dim=0; dim<pt->ndim; dim++)
            *(acceleration++) = 0;
        pt++;
    }
}

void clearAccelerations(chunk_particles_t* particlesChunk) {
    clearAccelerations(particlesChunk->particles, particlesChunk->nParticle);
}

void update_particles(particle_t* particles, int nParticle, double timeStep) {
    clearAccelerations(particles, nParticle);
    particle_t* one = particles;
    for(int i=0; i<nParticle; i++) {
        particle_t* another = particles;
        for(int j=0; j<nParticle; j++) {
            if (i!=j) 
                updateAcceleration(one, another);
            another++;
        }
        one++;
    }

    one = particles;
    for(int i=0; i<nParticle; i++) {
        updateVelocity(one, timeStep);
        updatePosition(one, timeStep);
    }
}

void update_particles(chunk_particles_t* particlesChunk, double timeStep) {
    update_particles(particlesChunk->particles, particlesChunk->nParticle, timeStep);
}
