#include "user_def/gravity_particle.h"


void gravityUpdateAcceleration(particle_t* one, particle_t* another) {
    const double G = 6.67e-11;
    double* mass1 = one->features;
    double* mass2 = another->features;

    double* positionDiff = nullptr;
    // printf("%d ", one->ndim);
    get_distance_vec(one->position, another->position, one->ndim, &positionDiff);
    double distant = get_length(positionDiff, one->ndim);
    distant = distant * distant * distant;
    double coeff1 = G / distant * (*mass2);
    double coeff2 = G / distant * (*mass1);
    double* acc = one->acceleration;
    double* positionDiff_ = positionDiff;
    for(int n=0; n<one->ndim; n++)
        *(acc++) += *(positionDiff_++) * coeff1;
    acc = another->acceleration;
    positionDiff_ = positionDiff;
    for(int n=0; n<another->ndim; n++)
        *(acc++) -= *(positionDiff_++) * coeff2;
    free(positionDiff);
}


