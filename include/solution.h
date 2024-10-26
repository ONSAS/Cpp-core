#ifndef SOLUTION
#define SOLUTION

#include <armadillo>
using namespace arma;

class Solution{
    private:
        double time;
    public:
        Solution();
        Solution(double time, vec disps);
        // Racional(const Solution &s); // constructor por copia
        // Racional suma(Racional &r);
        // bool igualdad(Racional &r);
        double get_time();
        vec get_disps();
        // friend bool operator ==(const Racional &r1, const Racional &r2);
        // friend Racional operator +(const Racional &r1, const Racional &r2);
        // operator == quiere decir que voy a sobrecargar el "=="
        ~Solution();
};
#endif