#include <vector>
#include <string>
#include <tuple>
#include <map>
#include <set>

#include "eigen-3.4.0/Eigen/Dense"

#include "basis.h"
#include "helper.h"

#pragma once

struct Atom {
    int id;
    float z;
    Eigen::Vector3d coord;
    std::vector<BasisFunction *> basis;

    ~Atom();
};

class Molecule {
    public:
        void print_geom();
        void print_zeta();
        
        // pass pointers to geometry for optimizer
        std::vector<double *> get_geom();

        // number of electrons
        int nel, nel_alp, nel_bet;
        
        // stores atoms with all data
        std::vector<Atom *> atoms;

        Molecule(std::string const &path);
        ~Molecule();
};