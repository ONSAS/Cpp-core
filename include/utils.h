#ifndef UTILS
#define UTILS

#include <armadillo>
using namespace arma;

class Utils{
  public:
    ivec nodes2dofs( ivec nodes, int degreesPerNode );
}

#endif