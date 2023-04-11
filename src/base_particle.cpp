#include <vector>
#include <cstdint>
#include <sstream>
#include <iostream>
#include "base_particle.h"
using namespace std;

vector<double> base_particle::compute_acceleration(base_particle other) {

    cerr << "The Function is not implemented.";
    vector<double> a(1, 0);
    return a;
}

void base_particle::compute_acceleration(Container particles) {
    for (auto iter: particles.particles) {
        this->compute_acceleration(iter);
    }
}

void base_particle::update_velocity(double timestep) {
    for (int i = 0; i < velocity.size(); i++) {
        velocity[i] += timestep * acceleration[i];
    }
}

void base_particle::update_position(double timestep) {
    for (int i = 0; i < velocity.size(); i++) {
        position[i] += timestep * velocity[i];
    }
}