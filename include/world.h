#ifndef WORLD
#define WORLD

#include "base_particle.h"
#include <cstdint>
#include <sstream>
using namespace std;

class world {
    public:
        world() {}
    private:
        double L0, L1, W0, W1, H0, H1;
        vector<base_particle> all_particle;
};

#endif