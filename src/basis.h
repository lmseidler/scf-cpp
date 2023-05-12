#include <vector>

#include "eigen-3.4.0/Eigen/Dense"

#pragma once

class BasisFunction {
    public:
        // id of the atom to which the basis function belongs to
        int atid;

        // id of the basis function in the set of all basis function
        int id;

        // quantum numbers for basis function
        int n;
        int l;

        // store the zeta value for optimization
        double zeta;

        // coefficients, normalization constants and exponents of primitives
        std::vector<double> d_set;
        std::vector<double> a_set;
        
        void expand_slater(int nprim);

        BasisFunction(int atid, int n, int l, double zeta);
};