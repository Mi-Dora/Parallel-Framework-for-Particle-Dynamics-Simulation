//
// Created by 杜煜 on 2023/4/10.
//
# include "utils.h"

vector<double> get_distance_vec(vector<double> &original, vector<double> &target){
    vector<double> dis_vec(3, 0);
    dis_vec[0] = target[0] - original[0];
    dis_vec[1] = target[1] - original[1];
    dis_vec[2] = target[2] - original[2];
    return dis_vec;
}

inline double get_length(vector<double> &vec){
    return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

void l2_norm(vector<double> &vec, double &vec_length){
    vec[0] /= vec_length;
    vec[1] /= vec_length;
    vec[2] /= vec_length;
}

