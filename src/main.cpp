
#include <iostream>
#include <armadillo>

#include "solution.h"

using namespace std  ;
using namespace arma ;

int main(){

    // read model
    double final_time = 2.0, delta_time = .5;
    printf(" %12.2e ", final_time);
    printf(" %12.2e ", delta_time);

    // vec U(2); U.ones();
    // Solution sol = Solution(0.0, U);

    // cout << " U " << sol.get_disps() << endl;

    // // compute time zero solution
    
    // // solve times solutions
    // double curr_time = 0.0, next_time = delta_time ;
    // while ( next_time < final_time  ){
    //     printf(" %12.2e ", next_time);
    //     curr_time = next_time ;
    //     next_time = curr_time + delta_time ; 
    // }
    // // generate visualization

    return 0;
}