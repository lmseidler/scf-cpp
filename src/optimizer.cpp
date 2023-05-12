#include <iostream>
#include <stdexcept>
#include <chrono>

#include "optimizer.h"

// parameter to adjust for steepest descent
// TODO - could still be tweaked
const double ETA = 1.0;
const double THRESHOLD_NORM = 0.003;
const double STEPSIZE_GEOM = 0.1;
const double STEPSIZE_BAS = 0.1;
const int MAX_STEPS = 100;

// #######################################
// Optimizer
// #######################################

// ######
// PUBLIC
// ######

bool Optimizer::run() {
    return steepest_descent();
}

Optimizer::Optimizer(Molecule &mol, int solver_type, double solver_threshold) : 
rhf(mol, solver_threshold), 
uhf(mol, solver_threshold), 
solver(nullptr), 
mol(mol) {
    // set solver type dyamically
    switch (solver_type) {
        case 0:
            solver = &rhf;
            break;
        case 1:
            solver = &uhf;
            break;
        default:
            solver = nullptr;
            break;
    }
}

Optimizer::~Optimizer() {}

// #########
// PROTECTED
// #########

// optimizer objective function definition
// should be defined for every optimizer separately
double Optimizer::objective() {
    return 0.0;
}

// function to print current state of optimization
void Optimizer::callback() {}

// function to print out final results
void Optimizer::results() {}

// optimize parameter set with steepest descent algorithm
bool Optimizer::steepest_descent() {
    // calculate initial guess gradient with arbitrary delta
    double stepsize_orig = stepsize;
    grad();
    std::cout << "Gradient norm: " << g.norm() << std::endl;

    int steps = 0;
    // while gradient norm above threshold and less than max allowed steps taken:
    while (g.norm() > THRESHOLD_NORM && steps < MAX_STEPS) {
        // update parameter set
        for (int i = 0; i < para.size(); i++) {
            *para[i] += ETA * g(i);
        }
        
        // calculate new gradient
        grad();
        
        // recalculate gradient if norm is too high (leading to too large step)
        for (int tries = 0; g.norm() > 5.0 && tries < MAX_STEPS; tries++) {
            stepsize *= 0.5;
            grad();
        }

        // reset stepsize
        stepsize = stepsize_orig;

        // run callback function to print results etc.
        callback();

        std::cout << "Gradient norm: " << g.norm() << "\n" << std::endl;
        steps++;
    }
    
    // final report for the optimization
    results();

    // return success of optimization
    if (steps == MAX_STEPS) {
        return false;
    } else {
        return true;
    }
}

// calculates the gradient by making two objective function evaluations
void Optimizer::grad() {
    std::cout << "Calculating gradient..." << std::endl;
    // local variables for objective function evaluations for every parameter
    double c1, c2;
    for (int i = 0; i < para.size(); i++) {
        // increment parameter
        *para[i] += stepsize;
        c1 = objective();

        // decrement parameter
        *para[i] -= 2 * stepsize;
        c2 = objective();

        // back to original value
        *para[i] += stepsize;

        // calculate gradient component for parameter
        g(i) = (c2 - c1) / (2 * stepsize);
    }
}

// #######################################
// GeomOptimizer
// #######################################

// ######
// PUBLIC
// ######

GeomOptimizer::GeomOptimizer(Molecule &mol, int solver_type, double solver_threshold) : Optimizer(mol, solver_type, solver_threshold) {
    // get parameters from mol
    this->para = mol.get_geom();

    // initialize gradient
    g = Eigen::VectorXd(para.size());

    // set stepsize
    stepsize = STEPSIZE_GEOM;

    // print initial energy
    std::cout << "Initial energy: " << objective() << std::endl;
}

GeomOptimizer::~GeomOptimizer() {}

// #########
// PROTECTED
// #########

// evaluate objective function
double GeomOptimizer::objective() {
    // run solver with updated molecule
    if (solver != nullptr) {
        return solver->run(false);
    } else {
        throw(std::invalid_argument("Solver for GeomOptimizer not initialized correctly."));
    }
}

void GeomOptimizer::callback() {
    std::cout << "Current geometry: " << std::endl;
    mol.print_geom();
    std::cout << "Current energy: " << objective() << std::endl;
}

void GeomOptimizer::results() {
    std::cout << "Final geometry: " << std::endl;
    mol.print_geom();
    std::cout << "Final energy: " << objective() << std::endl;
}

// #######################################
// BasisOptimizer
// #######################################

// ######
// PUBLIC
// ######

BasisOptimizer::BasisOptimizer(Molecule &mol, int solver_type, double solver_threshold) : Optimizer(mol, solver_type, solver_threshold) {
    // get parameters from molecule
    this->para = solver->get_basis();

    // initialize gradient
    g = Eigen::VectorXd(para.size());

    // set stepsize
    stepsize = STEPSIZE_BAS;

    // print initial zeta values
    mol.print_zeta();

    // print initial energy
    std::cout << "Initial energy: " << objective() << std::endl;
}

BasisOptimizer::~BasisOptimizer() {}

// #########
// PROTECTED
// #########

// evaluate objective function
double BasisOptimizer::objective() {
    // run solver with updated basis
    if (solver != nullptr) {
        // reinitialize basis functions
        solver->reinit();

        return solver->run(false);
    } else {
        throw(std::invalid_argument("Solver for BasisOptimizer not initialized correctly."));
    }
}

void BasisOptimizer::callback() {
    std::cout << "Current zeta values: " << std::endl;
    mol.print_zeta();
    std::cout << "Current energy: " << objective() << std::endl;
}

void BasisOptimizer::results() {
    std::cout << "Final zeta values: " << std::endl;
    mol.print_zeta();
    std::cout << "Final energy: " << objective() << std::endl;
}