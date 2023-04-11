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

vector<double> get_distance_vec(vector<double> &original, vector<double> &target);

double get_length(vector<double> &vec);

void l2_norm(vector<double> &vec, double &length);

#endif //PARALLEL_FRAMEWORK_FOR_PARTICLE_DYNAMICS_SIMULATION_UTILS_H
