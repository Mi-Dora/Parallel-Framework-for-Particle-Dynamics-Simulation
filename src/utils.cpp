//
// Created by 杜煜 on 2023/4/10.
//
# include "utils.h"
#include "stdlib.h"

// vector<double> get_distance_vec(vector<double> &original, vector<double> &target){
    // vector<double> dis_vec(3, 0);
    // dis_vec[0] = target[0] - original[0];
    // dis_vec[1] = target[1] - original[1];
    // dis_vec[2] = target[2] - original[2];
    // return dis_vec;
// }

inline void get_distance_vec(const double* vecA, const double* vecB, const int ndim, double** vecRet) {
    double* ret;
    double* a = const_cast<double*>(vecA);
    double* b = const_cast<double*>(vecB);
    ret = static_cast<double*>(malloc(sizeof(double)*ndim));
    for(size_t n=0; n<ndim; n++) 
        *(ret++) = *(a++) - *(b++);
    *vecRet = ret;
}


inline double get_length(const double* vec, const int ndim) {
    double ret = 0.0;
    double* _vec = const_cast<double*>(vec);
    for(size_t n=0; n<ndim; n++) {
        ret += (*_vec) * (*_vec++);
    }
    return sqrt(ret);
}

inline void l2_norm(const double* vec, const int ndim, double** vecRet) {
    double* ret;
    ret = static_cast<double*>(malloc(sizeof(double)*ndim));
    double length = get_length(vec, ndim);
    double* vec_ = const_cast<double*>(vec);
    for(size_t n=0; n<ndim; n++)
        *(ret++) = *(vec_) / length;
    *vecRet = ret;
}

inline void vec_add(const double* a, const double* b, const int ndim, double** sum) {
    double* a_ = const_cast<double*>(a);
    double* b_ = const_cast<double*>(b);
    double* sum_ = static_cast<double*>(malloc(sizeof(double)*ndim));

    for(size_t n=0; n<ndim; n++) 
        *(sum_++) = *(a_++) + *(b_++);
    *sum = sum_;
}

inline void vec_sub(const double* a, const double* b, const int ndim, double** ret) {
    double* a_ = const_cast<double*>(a);
    double* b_ = const_cast<double*>(b);
    double* ret_ = static_cast<double*>(malloc(sizeof(double)*ndim));

    for(size_t n=0; n<ndim; n++)
        *(ret_++) = *(a_++) - *(b_++);
    *ret = ret_;
}

