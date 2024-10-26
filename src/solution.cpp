
#include "solution.h"
#include <armadillo>

using namespace arma;

Solution::Solution(){
};

Solution::Solution(double t, vec U){
    float time = t;
    vec disps=U;
};

double Solution::get_time(){
    return time;
};

vec Solution::get_disps(){
    return disps;
};

Solution::~Solution(){
};
