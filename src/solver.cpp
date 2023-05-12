#include <iostream>
#include <tuple>
#include <cmath>

#include "eigen-3.4.0/Eigen/Dense"

#include "solver.h"
#include "molecule.h"

const int MAX_CYC = 50;

long TWOINT_TIMINGS = 0.0;
long NEWFOCK_TIMINGS = 0.0;

// #######################################
// Solver
// #######################################

// ######
// PUBLIC
// ######

double Solver::run(bool verbose) {
    return 0.0;
}

// reinitialize all basis functions
void Solver::reinit() {
    for (auto &&bf : basis) {
        bf->expand_slater(6);
    }
}

// collect pointers to basis function zeta values to be passed to optimizer
// using bf->id to ensure correct mapping
std::vector<double *> Solver::get_basis() {
    std::vector<double *> para = std::vector<double *>(nbas);

    for (const auto &bf : basis) {
        para[bf->id] = &(bf->zeta);
    }

    return para;
}

// collects basis functions from molecule and also calculates 1e and 2e integrals
Solver::Solver(const Molecule &mol, double threshold) : mol(mol) {
    this->threshold = threshold;

    // collect pointers to all basis functions in one vector
    this->basis = std::vector<BasisFunction *>();
    nbas = 0;
    for (Atom *atom : mol.atoms) {
        for (BasisFunction *bf : atom->basis) {
            basis.push_back(bf);
            bf->id = nbas;
            nbas++;
        }
    }
}

Solver::~Solver() {}

// #########
// PROTECTED
// #########

// calculates nuclear-nuclear interactions
void Solver::nuclear() {

    // total interaction energy
    v_nn = 0.0;

    // loop over all unique atom pairs
    for (int i = 0; i < mol.atoms.size()-1; i++) {
        for (int j = i+1; j < mol.atoms.size(); j++) {
            v_nn += mol.atoms[i]->z * mol.atoms[j]->z / std::sqrt(distance(mol.atoms[i]->coord, mol.atoms[j]->coord));
        }
    }
}

// executes the loop over the primitives of two CGTOs a and b
// returns the overlap, kinetic energy and nuclear attraction
std::array<double, 3> Solver::oneint_prim(int a, int b) {

    // local variables for readability
    Eigen::Vector3d auf_a = mol.atoms[basis[a]->atid]->coord;
    Eigen::Vector3d auf_b = mol.atoms[basis[b]->atid]->coord;
    double r_ab = distance(auf_a, auf_b);
    double s = 0;
    
    std::array<double, 3> stv = {0.0, 0.0, 0.0};
    
    // loop over all primitives
    for (int i = 0; i < basis[a]->a_set.size(); i++) {
        for (int j = 0; j < basis[b]->a_set.size(); j++) {

            // local variables for readability
            double a_i = basis[a]->a_set[i];
            double a_j = basis[b]->a_set[j];
            double a_p = a_i + a_j; 
            double d_ij =  basis[a]->d_set[i] * basis[b]->d_set[j];

            // calculate overlap
            s = d_ij * std::exp(-r_ab * (a_i * a_j) / a_p) * std::pow(std::sqrt(M_PI / a_p), 3);
            stv[0] += s;

            // calculate kinetic energy
            stv[1] += (a_i * a_j) / a_p * (3.0 - 2.0 * r_ab * (a_i * a_j) / a_p) * s;

            // calculate nuclear attraction
            // loop over all atoms
            // mixed gaussian at new aufpunkt
            Eigen::Vector3d auf_p = (a_i * auf_a + a_j * auf_b) / a_p;
            for (const Atom *atom : mol.atoms) {
                stv[2] -=   d_ij * 2 * M_PI / a_p 
                            * std::exp(-r_ab * (a_i * a_j) / a_p) 
                            * atom->z * boysf0(
                                a_p * distance(auf_p, atom->coord)
                            );
            }
        }
    }

    return stv;
}

// computes hcore and overlap matrix
// secondary objective for optimizations
void Solver::oneint() {
    // setup half-vectors for hcore and overlap matrix
    Eigen::VectorXd hcore_compact = Eigen::VectorXd::Zero(binom(nbas + 1, 2));
    Eigen::VectorXd s_compact = Eigen::VectorXd::Zero(binom(nbas + 1, 2));
    
    // VERSION1: 
    // loop over all atoms
        // loop over all unique basis set function pairs within each atom
            // calculate overlap and kinetic energy
                // loop over all atoms
                // calculate nuclear interaction

    // loop over all unique atom pairs
        // loop over all unique basis function pairs within the atom pair
            // calculate overlap and kinetic energy
            // loop over all atoms
                // calculate nuclear interaction
    // TODO - combine into one atom-loop
    // calculate h_core on the way, return it

    // VERSION2:
    // half-vector index (no mapping needed since later you can use hcore and s)
    int ab = 0;

    // loop over all unique basis function pairs within the molecule
    for (int a = 0; a < nbas; a++) {
        for (int b = a; b < nbas; b++) {
            // calculate overlap and kinetic energy and nuclear attraction
            std::array<double, 3> stv = oneint_prim(a, b);
            s_compact(ab) = stv[0];
            hcore_compact(ab) = stv[1] + stv[2];
            ab++;
        }
    }
    // unpack into hcore and overlap matrix
    hcore = unpack(hcore_compact);
    s = unpack(s_compact);
}

// executes the loop over all primitives of four CGTOs a, b, c, d
// returns the interaction energy
// TODO - implement vectorization (__m256)
// TODO - compiler flags
double Solver::twoint_prim(int a, int b, int c, int d) {

    // local variables for readability
    Eigen::Vector3d auf_a, auf_b, auf_c, auf_d;
    auf_a = mol.atoms[basis[a]->atid]->coord;
    auf_b = mol.atoms[basis[b]->atid]->coord;
    auf_c = mol.atoms[basis[c]->atid]->coord;
    auf_d = mol.atoms[basis[d]->atid]->coord;
    double r_ab = distance(auf_a, auf_b);
    double r_cd = distance(auf_c, auf_d);

    // integral value for (ab|cd)
    double g = 0.0;

    // loop over all primitives
    for (int i = 0; i < basis[a]->a_set.size(); i++) {
        for (int j = 0; j < basis[b]->a_set.size(); j++) {
            // local variables for readability
            double a_i = basis[a]->a_set[i];
            double a_j = basis[b]->a_set[j];
            double a_p = a_i + a_j; 
            double d_ij =  basis[a]->d_set[i] * basis[b]->d_set[j];
                
            // mixed gaussian at new aufpunkt
            Eigen::Vector3d auf_p = (a_i * auf_a + a_j * auf_b) / a_p;

            for (int k = 0; k < basis[c]->a_set.size(); k++) {
                for (int l = 0; l < basis[d]->a_set.size(); l++) {
                    // local variables for readability
                    double a_k = basis[c]->a_set[k];
                    double a_l = basis[d]->a_set[l];
                    double a_q = a_k + a_l; 
                    double d_kl =  basis[c]->d_set[k] * basis[d]->d_set[l];

                    // mixed gaussian at new aufpunkt
                    Eigen::Vector3d auf_q = (a_k * auf_c + a_l * auf_d) / a_q;

                    double r_pq = distance(auf_p, auf_q);

                    // evaluate integral
                    g +=    d_ij * d_kl 
                            * exp(-(r_ab * (a_i * a_j) / a_p) - (r_cd * (a_k * a_l) / a_q)) 
                            * 2 * std::pow(M_PI, 2.5) / (a_p * a_q * std::sqrt(a_p + a_q)) 
                            * boysf0(r_pq * (a_p * a_q) / (a_p + a_q));
                }
            }
        }
    }

    return g;
}

// calculates 2e integrals
// primary objective for optimizations
void Solver::twoint() {

    // setup index pairs by constructing all unique pairwise combinations of indices
    int nd = binom(nbas+1, 2);
    std::vector<std::pair<int, int>> pairs(nd);
    int ij = 0;
    for (int i = 0; i < nbas; i++) {
        for (int j = i; j < nbas; j++) {
            pairs[ij] = {i, j};
            ij++;
        }
    }

    // setup vector for 2e integral tensor
    g_abcd = Eigen::VectorXd::Zero(binom(nd + 1, 2));

    // index for vector: abcd
    int a, b, c, d, abcd = 0;

    // construct all unique combinations of pairs
    for (int i = 0; i < nd; i++) {
        for (int j = i; j < nd; j++) {
            // get basis function indices a, b, c, d from combinations of index pairs
            a = pairs[i].first;
            b = pairs[i].second;
            c = pairs[j].first;
            d = pairs[j].second;

            // calculate 2e integral (ab|cd)
            g_abcd(abcd) = twoint_prim(a, b, c, d);
            abcd_map[{{a, b}, {c, d}}] = abcd;
            abcd++;
        }
    }
}

// #######################################
// RHFSolver
// #######################################

// ######
// PUBLIC
// ######

RHFSolver::RHFSolver(const Molecule &mol, double threshold) : Solver(mol, threshold) {}

// runs RHF SCF procedure and returns final energy
double RHFSolver::run(bool verbose) {
    // calculate 1e and 2e integrals
    oneint();
    twoint();

    // calculate classical nuclear-nuclear interaction
    if (verbose) {
        std::cout << "\nCalculating nuclear-nuclear interaction..." << std::endl;
    }
    nuclear();
    if (verbose) {
        std::cout << "V_NN = " << v_nn << "\n" << std::endl;
    }

    // calculate electronic energy (initial guess)
    if (verbose) {
        std::cout << "Building initial guess..." << std::endl;
    }
    guess();
    if (verbose) {
        std::cout << "Initial guess density matrix:\n";

        for (int i = 0; i < p.rows(); i++) {
            for (int j = 0; j < p.cols(); j++) {
                printf("%-10.8f ", p(i, j));
            }
            std::cout << "\n";
        }

        std::cout << "\n";
    }
    double e = energy();

    e += v_nn;

    int n_cyc = 0;
    double e_prev, rmsd, delta_e;

    if (verbose) {
        printf("%-5s%-15s%-15s%-15s\n", "NCYC", "ETOT", "EL", "DELTA");
        printf("%-5d%-15.8f%-15.8f%-15.8f\n", n_cyc, e, e - v_nn, delta_e); 
    }
    

    do {
        n_cyc++;
        e_prev = e;

        // build new fock and orthogonalize
        // also build new density matrix
        new_fock();
        
        // get new energy
        e = energy();

        e += v_nn;
        delta_e = e - e_prev;

        if (verbose) {
            printf("%-5d%-15.8f%-15.8f%-15.8f\n", n_cyc, e, e - v_nn, delta_e);
        }

    } while (abs(delta_e) > threshold and n_cyc < MAX_CYC);

    if (verbose) {
        std::cout << "\nFinal SCF energy: " << e << "\n" << std::endl;

        std::cout << "Mulliken charges:" << std::endl;
        mulliken();
    }
    return e;
}

// #######
// PRIVATE
// #######

// calculate an initial density matrix from the ao coefficient matrix
void RHFSolver::dens_init() {
    // setup matrix with occupation numbers
    nocc = Eigen::MatrixXd::Zero(nbas, nbas);
    for (int i = 0; i < mol.nel / 2; i++) {
        nocc(i, i) = 2;
    }

    // setup density matrix
    p = Eigen::MatrixXd::Zero(nbas, nbas);
    p = ao_coeff * nocc * ao_coeff.transpose();
}


// provides the molecule with an initial guess for the density matrix,
// as well as fock, fock_ortho, overlap, ao_coeff, coeff, s_diag, s_evecs, ortho, fock_ortho
void RHFSolver::guess() {

    // build orthogonalization matrix
    std::tie(s_diag, s_evecs) = diag_nxn(s);
    ortho = s_evecs * s_diag.array().sqrt().matrix().inverse() * s_evecs.transpose();

    // initial guess for fock
    fock = hcore;
    fock_ortho = ortho.transpose() * fock * ortho;

    // get mo coefficient matrix from the orthogonalized fock matrix
    std::tie(mo_energies, coeff) = diag_nxn(fock_ortho);

    // get ao_coeff by orthogonalizing coeff
    ao_coeff = ortho * coeff;

    // setup guess for density matrix
    dens_init();
}

// evaluates the Fock equation for the electronic energy
double RHFSolver::energy() {
    return 0.5 * ((hcore + fock) * p).trace();
}

// build new density matrix from ao coefficients
void RHFSolver::dens() {
    p = ao_coeff * nocc * ao_coeff.transpose();
}

// creates new fock matrix including 2e integrals and orthogonalizes
void RHFSolver::new_fock() {
    // always start from hcore
    fock = hcore;

    // add 2e integrals multiplied by p
    for (int a = 0; a < nbas; a++) {
        for (int b = 0; b < nbas; b++) {
            for (int c = 0; c < nbas; c++) {
                for (int d = 0; d < nbas; d++) {
                    fock(a, b) += p(c, d) * (
                        g_abcd(abcd_map[{{a, b}, {c, d}}]) - 0.5 * g_abcd(abcd_map[{{a, c}, {b, d}}])
                    );
                }
            }
        }
    }

    // orthogonalize
    fock_ortho = ortho.transpose() * fock * ortho;

    // diagonalize
    std::tie(mo_energies, coeff) = diag_nxn(fock_ortho);

    // get ao_coeff by orthogonalizing coeff
    ao_coeff = ortho * coeff;

    // build new density matrix from new coefficients
    dens();
}

// perform a Mulliken population analysis
void RHFSolver::mulliken() {
    double qi;
    Eigen::MatrixXd p_s = p * s;

    for (int i = 0; i < mol.atoms.size(); i++) {
        qi = mol.atoms[i]->z;
        for (BasisFunction *bf : mol.atoms[i]->basis) {
            qi -= p_s(bf->id, bf->id);
        }
        std::cout << "Atom " << i << ": " << qi << std::endl;
    }

    std::cout << "\n" << std::endl;
}

// #######################################
// UHFSolver
// #######################################

// ######
// PUBLIC
// ######

UHFSolver::UHFSolver(const Molecule &mol, double threshold) : Solver(mol, threshold) {}

// runs RHF SCF procedure and returns final energy
double UHFSolver::run(bool verbose) {
    // calculate 1e and 2e integrals
    oneint();
    twoint();

    // calculate classical nuclear-nuclear interaction
    if (verbose) {
        std::cout << "\nCalculating nuclear-nuclear interaction..." << std::endl;
    }
    nuclear();
    if (verbose) {
        std::cout << "V_NN = " << v_nn << "\n" << std::endl;
    }

    // calculate electronic energy (initial guess)
    if (verbose) {
        std::cout << "Building initial guess..." << std::endl;
    }
    guess();
    if (verbose) {
        std::cout << "Initial guess density matrix (alpha):\n";

        for (int i = 0; i < p_alp.rows(); i++) {
            for (int j = 0; j < p_alp.cols(); j++) {
                printf("%10.8f ", p_alp(i, j));
            }
            std::cout << "\n";
        }

        std::cout << "\n";

        std::cout << "Initial guess density matrix (beta):\n";

        for (int i = 0; i < p_bet.rows(); i++) {
            for (int j = 0; j < p_bet.cols(); j++) {
                printf("%-10.8f ", p_bet(i, j));
            }
            std::cout << "\n";
        }

        std::cout << "\n";
    }
    double e = energy();

    e += v_nn;

    int n_cyc = 0;
    double e_prev, rmsd, delta_e;

    if (verbose) {
        printf("%-5s%-15s%-15s%-15s\n", "NCYC", "ETOT", "EL", "DELTA");
        printf("%-5d%-15.8f%-15.8f%-15.8f\n", n_cyc, e, e - v_nn, delta_e); 
    }
    

    do {
        n_cyc++;
        e_prev = e;

        // build new fock and orthogonalize
        // also build new density matrix
        new_fock();
        
        // get new energy
        e = energy();

        e += v_nn;
        delta_e = e - e_prev;

        if (verbose) {
            printf("%-5d%-15.8f%-15.8f%-15.8f\n", n_cyc, e, e - v_nn, delta_e);
        }

    } while (abs(delta_e) > threshold and n_cyc < MAX_CYC);

    if (verbose) {
        std::cout << "\nFinal SCF energy: " << e << "\n" << std::endl;

        // calculate and print out spin contamination
        double spin_exp = (mol.nel_alp - mol.nel_bet) / 2.0 * ((mol.nel_alp - mol.nel_bet) / 2.0 + 1.0);
        double spin_calc = spin_exp + contamination();

        std::cout << "\nSpin expectation value: " << spin_exp << std::endl;
        std::cout << "Calculated spin: " << spin_calc << std::endl;
    }
    return e;
}

// #######
// PRIVATE
// #######

// creates new fock matrix including 2e integrals and orthogonalizes
void UHFSolver::new_fock() {
    // always start from hcore
    fock_alp = hcore;
    fock_bet = hcore;

    // add 2e integrals multiplied by p
    for (int a = 0; a < nbas; a++) {
        for (int b = 0; b < nbas; b++) {
            for (int c = 0; c < nbas; c++) {
                for (int d = 0; d < nbas; d++) {
                    // fock matrix for alpha spin
                    fock_alp(a, b) += (p_alp(c, d) + p_bet(c, d)) * g_abcd(abcd_map[{{a, b}, {c, d}}]);
                    fock_alp(a, b) -= p_alp(c, d) * g_abcd(abcd_map[{{a, c}, {b, d}}]);

                    // fock matrix for beta spin
                    fock_bet(a, b) += (p_alp(c, d) + p_bet(c, d)) * g_abcd(abcd_map[{{a, b}, {c, d}}]);
                    fock_bet(a, b) -= p_bet(c, d) * g_abcd(abcd_map[{{a, c}, {b, d}}]);
                }
            }
        }
    }

    // orthogonalize
    fock_alp_ortho = ortho.transpose() * fock_alp * ortho;
    fock_bet_ortho = ortho.transpose() * fock_bet * ortho;

    // diagonalize
    std::tie(mo_energies_alp, coeff_alp) = diag_nxn(fock_alp_ortho);
    std::tie(mo_energies_bet, coeff_bet) = diag_nxn(fock_bet_ortho);

    // get ao_coeff by orthogonalizing coeff
    ao_coeff_alp = ortho * coeff_alp;
    ao_coeff_bet = ortho * coeff_bet;

    // build new density matrices from new coefficients
    dens();
}

// build new density matrices from ao coefficients
void UHFSolver::dens() {
    p_alp = ao_coeff_alp * nocc_alp * ao_coeff_alp.transpose();
    p_bet = ao_coeff_bet * nocc_bet * ao_coeff_bet.transpose();
}

// evaluates the Fock equation for the electronic energy
double UHFSolver::energy() {
    double e = 0.0;
    for (int a = 0; a < nbas; a++) {
        for (int b = 0; b < nbas; b++) {
            e += 0.5 * (p_alp(a, b) * (hcore(a, b) + fock_alp(a, b)) + p_bet(a, b) * (hcore(a, b) + fock_bet(a, b)));
        }
    }

    return e;
}

// calculate an initial density matrix from the ao coefficient matrix
void UHFSolver::dens_init() {
    // setup matrix with occupation numbers
    nocc_alp = Eigen::MatrixXd::Zero(nbas, nbas);
    nocc_bet = Eigen::MatrixXd::Zero(nbas, nbas);

    // alpha occupations
    for (int i = 0; i < mol.nel_alp; i++) {
        nocc_alp(i, i) = 1;
    }

    // beta occupations
    for (int i = 0; i < mol.nel_bet; i++) {
        nocc_bet(i, i) = 1;
    }

    // setup density matrix
    p_alp = Eigen::MatrixXd::Zero(nbas, nbas);
    p_bet = Eigen::MatrixXd::Zero(nbas, nbas);

    p_alp = ao_coeff_alp * nocc_alp * ao_coeff_alp.transpose();
    p_bet = ao_coeff_bet * nocc_bet * ao_coeff_bet.transpose();
}


// provides the molecule with an initial guess for the density matrix,
// as well as fock, fock_ortho, overlap, ao_coeff, coeff, s_diag, s_evecs, ortho, fock_ortho
void UHFSolver::guess() {

    // build orthogonalization matrix
    std::tie(s_diag, s_evecs) = diag_nxn(s);
    ortho = s_evecs * s_diag.array().sqrt().matrix().inverse() * s_evecs.transpose();

    // initial guess for fock
    fock_alp = hcore;
    fock_bet = hcore;
    fock_alp_ortho = ortho.transpose() * fock_alp * ortho;
    fock_bet_ortho = ortho.transpose() * fock_bet * ortho;

    // diagonalize
    std::tie(mo_energies_alp, coeff_alp) = diag_nxn(fock_alp_ortho);
    std::tie(mo_energies_bet, coeff_bet) = diag_nxn(fock_bet_ortho);

    // get ao_coeff by orthogonalizing coeff
    ao_coeff_alp = ortho * coeff_alp;
    ao_coeff_bet = ortho * coeff_bet;


    // setup guess for density matrices
    dens_init();
}

// calculate spin contamination value
// TODO - there seems to be a slight deviation from reference values
double UHFSolver::contamination() {
    double cont = mol.nel_bet;

    for (int i = 0; i < mol.nel_alp; i++) {
        for (int j = 0; j < mol.nel_bet; j++) {
            cont -= std::pow(coeff_alp.col(i).dot(coeff_bet.col(j)), 2);
        }
    }

    return cont;
}