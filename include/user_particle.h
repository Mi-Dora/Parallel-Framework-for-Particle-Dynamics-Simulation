#ifndef USER_PARTICLE
#define USER_PARTICLE

#include <vector>
#include <cstdint>
#include <sstream>
#include "base_particle.h"
using namespace std;

class user_particle : public base_particle {
    public:
        user_particle(vector<double>& position, vector<double>& acceleration, vector<double>& velocity, double mass, double charge)
        : base_particle(position, acceleration, velocity), mass(mass), charge(charge){}
        vector<double> compute_acceleration(user_particle& other);

    private:
        double mass;
        double charge;
        static double ke;
};

#endif