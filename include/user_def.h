#ifndef USER_DEF_H
#define USER_DEF_H

#include <vector>
#include <cstdint>
#include <sstream>
using namespace std;

class Particle{
    public:
        Particle(double l, double w, double h, double charge, double mass);
        Particle(vector<double> position, vector<double> force,
                vector<double> speed, double charge, double mass);

        vector<double> force_compute(Particle other);
        void update_force(vector<double> force);
        void update_speed(double timestep);
        void update_position();
    private:
        vector<double> position;
        vector<double> force;
        vector<double> speed;
        double charge;
        double mass;
};

#endif