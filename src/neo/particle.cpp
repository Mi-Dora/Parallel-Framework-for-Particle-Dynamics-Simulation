#include "neo/particle.h"
#include <vector>
#include <cassert>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cmath>

void updateVelocity(particle_t* particle, double timeStep) {
    for(int n=0; n<PARTICLE_N_DIM; n++)
        particle->velocity[n] += particle->acceleration[n]*timeStep;
}

void updatePosition(particle_t* particle, double timeStep) {
    for(int n=0; n<PARTICLE_N_DIM; n++)
        particle->position[n] += particle->velocity[n]*timeStep;
}

void clearAcceleration(particle_t* particle) {
    for(int n=0; n<PARTICLE_N_DIM; n++)
        particle->acceleration[n] = 0.0;
}

void loadParticles(const std::string filename, std::vector<particle_t>* particles) {
    assert(particles->size() == 0);
    std::ifstream file(filename, std::ios::in);
    if(!file.is_open()) {
        std::fprintf(stderr, "[%s:%d]: can not open the file %s\n", __FILE__, __LINE__, filename);
        exit(EXIT_FAILURE);
    }

    int nParticles, nDim;
    std::string line;
    std::getline(file, line);
    std::sscanf(line.c_str(), "%d %d", &nParticles, &nDim);
    assert(nDim == PARTICLE_N_DIM);
    std::stringstream lineStream;
    int n = 0;
    char* bar;
    double foo;
    particles->resize(nParticles);
    while(std::getline(file, line) && n<nParticles) {
        particle_t& particle = particles->at(n);
        bar = strtok(const_cast<char*>(line.c_str()), " ");
        foo = atof(bar);
        particle.id = static_cast<uint64_t>(foo);
        for(int nn=0; nn<PARTICLE_N_DIM; nn++) {
            bar = strtok(NULL, " ");
            particle.position[nn] = atof(bar);
        }
        for(int nn=0; nn<PARTICLE_N_DIM; nn++) {
            bar = strtok(NULL, " ");
            particle.velocity[nn] = atof(bar);
        }
        for(int nn=0; nn<PARTICLE_N_DIM; nn++) {
            bar = strtok(NULL, " ");
            particle.acceleration[nn] = atof(bar);
        }
        bar = strtok(NULL, " ");
        foo = atof(bar);
        int nFeat = static_cast<int>(foo);
        assert(nFeat < PARTICLE_N_FEAT);
        for(int nn=0; nn<nFeat; nn++) {
            bar = strtok(NULL, " ");
            particle.features[nn] = atof(bar);
        }
        for(int nn=nFeat; nn<PARTICLE_N_FEAT; nn++) {
            particle.features[nn] = nan("");
        }
        particle.enabled=true;
        n++;
    }
    file.close();
    if(n != nParticles) {
        std::fprintf(stderr, "[%s:%d]: unmatched particle number %d vs %d\n", __FILE__, __LINE__, n, nParticles);
        exit(EXIT_FAILURE);
    }
}

void storeParticles(const std::string filename, std::vector<particle_t>* particles) {
    std::ofstream file(filename, std::ios::out|std::ios::trunc);
    if(!file.is_open()) {
        std::fprintf(stderr, "[%s:%d]: can not open the file %s\n", __FILE__, __LINE__, filename);
        exit(EXIT_FAILURE);
    }
    file << particles->size() << " " << PARTICLE_N_DIM << std::endl;
    file.flush();
    for(int n=0; n<particles->size(); n++) {
        particle_t& particle = particles->at(n);
        file << std::scientific << particle.id << " ";
        for(int nn=0; nn<PARTICLE_N_DIM; nn++)
            file << std::scientific << particle.position[nn] << " ";
        for(int nn=0; nn<PARTICLE_N_DIM; nn++)
            file << std::scientific << particle.velocity[nn] << " ";
        for(int nn=0; nn<PARTICLE_N_DIM; nn++)
            file << std::scientific << particle.acceleration[nn] << " ";
        for(int nn=0; nn<PARTICLE_N_FEAT; nn++) {
            file << std::scientific << particle.features[nn];
            if(nn+1 != PARTICLE_N_FEAT)
                file << " ";
            else
                file << std::endl;
        }
    }
    file.close();
}
