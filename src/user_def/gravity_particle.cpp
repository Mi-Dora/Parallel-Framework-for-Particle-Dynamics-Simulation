#include "user_def/gravity_particle.h"
#include "utils.h"


inline void gravityUpdateAcceleration(particle_t* one, particle_t* another) {
    const double G = 6.67e-11;
    double* mass1 = one->features;
    double* mass2 = another->features;

    double** positionDiff;
    get_distance_vec(one->position, another->position, one->ndim, positionDiff);
    double distant = get_length(*positionDiff, one->ndim);
    distant = distant * distant * distant;
    distant = G / distant * (*mass2);
    double* acc = one->acceleration;
    double* positionDiff_ = *positionDiff;
    for(size_t n=0; n<one->ndim; n++)
        *(acc++) += *(positionDiff_++) * distant;
    free(*positionDiff);
}


