#ifndef BASE_PARTICLE
#define BASE_PARTICLE

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
        base_particle(vector<double>& position, vector<double>& acceleration, vector<double>& velocity):
        position(position), acceleration(acceleration), velocity(velocity){}
        vector<double> compute_acceleration(base_particle other);
        void compute_acceleration(Container particles);
        void update_velocity(double timestep);
        void update_position(double timestep);
    protected:
        vector<double> position;
        vector<double> acceleration;
        vector<double> velocity;
};

#endif