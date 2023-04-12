//
// Created by Yu Du on 2023/4/10.
//

#ifndef PARALLEL_FRAMEWORK_FOR_PARTICLE_DYNAMICS_SIMULATION_UTILS_H
#define PARALLEL_FRAMEWORK_FOR_PARTICLE_DYNAMICS_SIMULATION_UTILS_H
#include <vector>
#include <cstdint>
#include <sstream>
#include <cmath>
using namespace std;

void get_distance_vec(const double* vecA, const double* vecB, const int ndim, double** vecRet);

double get_length(const double* vec, const int ndim);

void l2_norm(const double* vec, const int ndim, double** vecRet);

#endif //PARALLEL_FRAMEWORK_FOR_PARTICLE_DYNAMICS_SIMULATION_UTILS_H
