#include <vector>
#include <cstdint>
#include <sstream>
#include <user_def.h>
using namespace std;

Particle::Particle(double l, double w, double h, double charge, double mass): charge(charge), mass(mass) {
    position = {rand() * l, rand() * w, rand() * h};
    force = {0, 0, 0};
    speed = {0, 0, 0};
}

Particle::Particle(vector<double> position, vector<double> force,
                vector<double> speed, double charge, double mass): position(position),
                force(force), speed(speed), charge(charge), mass(mass) {}

void Particle::update_force(vector<double> force) {
    force = force;
}

void Particle::update_speed(double timestep) {
    speed[0] += force[0] * timestep / mass;
    speed[1] += force[1] * timestep / mass;
    speed[2] += force[2] * timestep / mass;

    // coefficient = timestep/mass;
    // speed += force * coefficient ?????
}

void Particle::update_position() {

}