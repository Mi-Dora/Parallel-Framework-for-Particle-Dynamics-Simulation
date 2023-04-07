#ifndef USER_DEF_H
#define USER_DEF_H

#include <vector>
#include <cstdint>
#include <sstream>
using namespace std;

class Container {
    public:
        vector<base_particle> particles;
};

class base_particle {
    public:
        vector<double> compute_acceleration(base_particle other);
        vector<double> compute_acceleration(Container particles);
        void update_velocity(double timestep);
        void update_position(double timestep);
    private:
        vector<double> position;
        vector<double> acceleration;
        vector<double> velocity;
};

#endif