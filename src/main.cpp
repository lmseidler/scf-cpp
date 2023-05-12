#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>

#include "eigen-3.4.0/Eigen/Dense"

#include "stringio.h"
#include "molecule.h"
#include "helper.h"
#include "solver.h"
#include "optimizer.h"

// solver energy change threshold
const double THRESHOLD_E = 0.0005;

// gets the setting for a given command line flag
char *get_cmd_option(char **begin, char **end, const std::string &flag) {
    char ** itr = std::find(begin, end, flag);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

// looks for a given command line flag
bool cmd_option_exists(char** begin, char** end, const std::string &flag) {
    return std::find(begin, end, flag) != end;
}

int main(int argc, char *argv[]) {
    // vars for command line arguments
    bool rhf, uhf, geom, bas = false;
    const char *path;

    // check argc
    if (argc > 5) {
        std::cout << "Too many arguments." << std::endl;
        std::cout << "Usage: scf [-i <path/to/input>] [-rhf / -uhf] [-geom / -bas]" << std::endl;
        return 1;
    } else {
        // parse command line arguments
        if (cmd_option_exists(argv, argv + argc, "-i")) {
            path = get_cmd_option(argv, argv + argc, "-i");
        } else {
            std::cout << "No input file given." << std::endl;
            return 1;
        }

        // determine solver type
        if (cmd_option_exists(argv, argv + argc, "-rhf")) {
            rhf = true;
            uhf = false;
        } else if (cmd_option_exists(argv, argv + argc, "-uhf")) {
            uhf = true;
            rhf = false;
        } else {
            // run RHF by default
            rhf = true;
            uhf = false;
        }

        // determine if an optimization should be executed
        if (cmd_option_exists(argv, argv + argc, "-geom")) {
            geom = true;
            bas = false;
        } else if (cmd_option_exists(argv, argv + argc, "-bas")) {
            bas = true;
            geom = false;
        } else {
            // by default run no optimization
            geom = false;
            bas = false;
        }
    }

    // initialize the molecule
    // reading geometry and basis set from input file
    Molecule mol(path);

    std::cout << "Read geometry:" << std::endl;
    mol.print_geom();

    // if optimization requested:
    if (geom || bas) {
        bool success = false;
        Optimizer *optimizer;

        // select optimizer type
        if (geom) {
            // check number of atoms
            if (mol.atoms.size() > 1) {
                optimizer = new GeomOptimizer(mol, rhf ? 0 : 1, THRESHOLD_E);
            } else {
                std::cout << "Only one atom found in molecule, geometry optimization not possible." << std::endl;
                return 1;
            }
        } else if (bas) {
            optimizer = new BasisOptimizer(mol, rhf ? 0 : 1, THRESHOLD_E);
        }

        success = optimizer->run();

        if (success) {
                std::cout << "Optimization successful!" << std::endl;
            } else {
                std::cout << "Optimization failed!" << std::endl;
        }

        delete optimizer;
    } else {
        // no optimization, just SCF
        // initialize solver with molecule, run verbose SCF
        if (rhf) {
            RHFSolver solver(mol, THRESHOLD_E);
            solver.run(true);
        } else if (uhf) {
            UHFSolver solver(mol, THRESHOLD_E);
            solver.run(true);
        }
    }

    return 0;
}