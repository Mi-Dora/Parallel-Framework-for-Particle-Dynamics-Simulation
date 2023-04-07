#include <vector>
#include <cstdint>
#include <sstream>
#include "base_particle.h"
using namespace std;

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