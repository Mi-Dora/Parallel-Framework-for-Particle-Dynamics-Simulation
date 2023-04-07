#include <vector>
#include <cstdint>
#include <sstream>
#include "base_particle.h"
using namespace std;

void base_particle::compute_acceleration(base_particle other) {

}

void base_particle::update_velocity(double timestep) {
    speed[0] += force[0] * timestep / mass;
    speed[1] += force[1] * timestep / mass;
    speed[2] += force[2] * timestep / mass;

    // coefficient = timestep/mass;
    // speed += force * coefficient ?????
}

void base_particle::update_position() {

}