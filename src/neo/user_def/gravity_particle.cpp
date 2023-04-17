#include "neo/user_def/gravity_particle.h"
#include "neo/particle.h"


void gravityUpdateAcceleration(particle_t* one, particle_t* another) {
    const double G = 6.67e-11;
    double mass1 = one->features[0];
    double mass2 = another->features[0];

    double* positionDiff = nullptr;
    // printf("%d ", one->ndim);
    get_distance_vec(one->position, another->position, PARTICLE_N_DIM, &positionDiff);
    double distant = get_length(positionDiff, PARTICLE_N_DIM);
    distant = distant * distant * distant;
    double coeff1 = G / distant * mass2;
    double coeff2 = G / distant * mass1;
    double* acc = one->acceleration;
    double* positionDiff_ = positionDiff;
    for(int n=0; n<PARTICLE_N_DIM; n++)
        *(acc++) += *(positionDiff_++) * coeff1;
    acc = another->acceleration;
    positionDiff_ = positionDiff;
    for(int n=0; n<PARTICLE_N_DIM; n++)
        *(acc++) -= *(positionDiff_++) * coeff2;
    free(positionDiff);
}


