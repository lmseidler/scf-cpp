#include <stdexcept>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <tuple>
#include <string>
#include <array>

#include "eigen-3.4.0/Eigen/Dense"

#include "eigen-3.4.0/Eigen/src/Core/Matrix.h"
#include "molecule.h"
#include "stringio.h"
#include "basis.h"
#include "helper.h"

// notes: 
// h_core: tensor of order 2 (== matrix)
// g_abcd: tensor of order 4

// ######
// PUBLIC
// ######

// prints out the geometry of the molecule (atomic numbers, coordinates)
void Molecule::print_geom() {
    for (const Atom *atom : atoms) {
        std::cout << atom->z << " ";
        for (int j = 0; j < 3; j++) {
            printf("%10.8f ", atom->coord[j]);
        }
        std::cout << std::endl;
    }
    std::cout << "\n";
}

// prints out the zeta values for the STO basis functions
void Molecule::print_zeta() {
    for (const auto &atom : atoms) {
        std::cout << "Atom " << atom->id << ":" << std::endl;
        for (const auto &bf : atom->basis) {
            std::cout << bf->zeta << std::endl;
        }
        std::cout << "\n" << std::endl;
    }
}

std::vector<double *> Molecule::get_geom() {
    std::vector<double *> para;

    for (const auto &atom : atoms) {
        for (auto &x : atom->coord) {
            para.push_back(&x);
        }
    }

    return para;
}

// molecule constructor reads data from input file and allocates atoms and basis functions
Molecule::Molecule(const std::string &path) {
    // try to open input file
    std::ifstream infile;
    infile.open(path);

    if (infile.is_open()) {

        // local variables
        std::string line;
        int nbasis;
        
        // read the first line
        std::getline(infile, line);
        std::vector<std::string> line_split = split(line);

        // set number of electrons
        nel = std::stoi(line_split[1]);

        // set number of alpha and beta electrons by just reading the difference (fourth number in first line of input file)
        int nel_delta = std::stoi(line_split[3]);
        if ((nel_delta % 2 == 0 && nel % 2 == 0) || (nel_delta % 2 == 1 && nel % 2 == 1)) {
            nel_alp = (nel + nel_delta) / 2;
            nel_bet = nel - nel_alp;
        } else {
            throw(std::invalid_argument("Multiplicity not possible with given number of electrons."));
        }

        while (std::getline(infile, line)) {
            // split line into substrings
            line_split = split(line);

            // assumes that a line containing 5 entries follows (as in the input examples)
            if (line_split.size() != 5) {
                throw(std::invalid_argument("Format error in input file."));
            }

            // allocate new atom
            Atom *new_atom = new Atom;
            new_atom->id = atoms.size();

            // set atom coordinates
            for (int i = 0; i < 3; i++) {
                new_atom->coord(i) = std::stod(line_split[i]);
            }

            // set charge
            new_atom->z = std::stof(line_split[3]);

            // consume all following basis function exponents and expand into slater (STO-6G hardcoded)
            nbasis = std::stoi(line_split[4]);
            for (int i = 0; i < nbasis; i++) {
                std::getline(infile, line);

                // allocate new basis function
                BasisFunction *new_cgto = new BasisFunction(atoms.size(), 1, 0, std::stod(line));

                // expand STO into CGTO
                new_cgto->expand_slater(6);

                // put new basis function into atom
                new_atom->basis.push_back(new_cgto);
            }
            // add new atom to molecule
            atoms.push_back(new_atom);
        }

        infile.close();
    } else {
        throw(std::invalid_argument("Could not open input file."));
    }
}

Molecule::~Molecule() {
    for (Atom *atom : atoms) {
        delete atom;
    }
}

Atom::~Atom() {
    for (BasisFunction *bf : basis) {
        delete bf;
    }
}