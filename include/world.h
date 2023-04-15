#ifndef WORLD
#define WORLD

#include "particle.h"
#include <cstdint>
#include <sstream>
using namespace std;

void generator();

class world {
    public:
        world(double L0, double L1, double W0, double W1, double H0, double H1): \
            L0(L0), L1(L1), W0(W0), W1(W1), H0(H0), H1(H1) {}
    private:
        double L0, L1, W0, W1, H0, H1;
        vector<user_particle> all_particle;
};

#endif